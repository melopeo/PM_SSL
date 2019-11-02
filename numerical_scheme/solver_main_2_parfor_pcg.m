function v = solver_main_2_parfor_pcg(vector, A_cell,p,lambda,scheme_str)
%     Solves the system (I + \lambda*Lp)v = vector
% INPUT: 
%      : vector    : (vector) column vector
%      : A_cell    : (cell array) cell array with Laplacian matrices
%      : p         : (scalar) power (everything assumes that p<0 and -p = integer)
%      : lambda    : (scalar)
%      : diagShift : (scalar) diagonal shift used to compute matrices of A_cell
%      : scheme_str: (string) indicates how to compute the solution 'v'
%                  : options: 'cmg'       : uses conditional multigrid Laplacian solver
%                           : 'quadrature': uses quadrature approximation
%                           : 'explicit'  : uses explicit computation with matlab
%      : maxiter   : (scalar) maximum number of iterations (applies only for 'cmg'
%                    and 'quadrature'

% restoredefaultpath
% addpath(genpath('functions'))
% addpath(genpath('PM'))
% addpath(genpath('CMG')) % Combinatorial Multigrid is a solver

warning('off', 'MATLAB:gmres:tooSmallTolerance');  
warning('off', 'MATLAB:pcg:tooSmallTolerance');  

numNodes = length(vector);
numLayers = length(A_cell);

switch scheme_str
    case 'pcg'
        % Ichol preconditioners
        threshold = 1e-4;
        L_array = ichol_data(A_cell, threshold);  
end

% % % % Define Mp times a vector function handle
switch scheme_str
    case 'pcg'
       f_Mp = @(x)  Mp_x_vector_pcg(A_cell,L_array,x,p);
    case 'explicit'
       Mp = get_Mp(A_cell,p);
       M  = eigs(Mp,1,'SM');
       m  = eigs(Mp,1,'LM');
end

% % % % Compute quadrature parameters
switch scheme_str
    case 'pcg'
        % % % % Estimate minimum and maximum eigenvalues of Mp
        M = eigs(f_Mp,numNodes,1,'LM');
        m = estimate_smallest_eigenvalue_of_Mp_weyls_theorem(A_cell,p);
        tol_contour  = 1.e-8;
        N            = 2*ceil(abs(log(tol_contour)*(log(M/m)+3)/(2*pi*pi)));
        [z,dzdt,K,k] = contour_points(m,M,N);
end

% % % % Main computation
switch scheme_str
    case 'pcg'
        gmres_maxiter = 5*abs(p);
        f_Lp  = @(x) Lp_x_vector(f_Mp,x,numLayers,p,m,M,N,z,dzdt,K,k,gmres_maxiter);
        v = solve_I_plus_sigmaLp(f_Lp,lambda,vector);
    case 'explicit'
        Mp = get_Mp(A_cell,p);
        Lp = length(A_cell)^(-1/p) * full(Mp)^(1/p); 
        v  = (eye(numNodes,numNodes) + lambda*Lp)\vector;
end


warning('on', 'MATLAB:gmres:tooSmallTolerance');  
warning('on', 'MATLAB:pcg:tooSmallTolerance'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function m = estimate_smallest_eigenvalue_of_Mp_weyls_theorem(A_array,p)
    m = 0;
    for i = 1 : length(A_array)
        m = m + eigs(A_array{i},1,'LM')^p;
    end

function Mp = get_Mp(A_array,p)
   Mp = zeros(size(A_array{1}));
   for i = 1 : length(A_array)
       Mp = Mp + full(A_array{i})^p;
   end
   
function v = solve_I_plus_sigmaLp(f_Lp, sigma, vector)
    f = @(x) x + sigma * f_Lp(x);
    
    tol = 1e-3; maxiter = 10;
    
    [v,~] = gmres(f, vector, [],  tol, maxiter);

function v = Lp_x_vector(f_Mp,vector,ell,p,m,M,N,z,dzdt,K,k,maxiter) 
    f = @(x) x^(1/p);
    vector = f_Mp(vector);
    v = zeros(size(vector));
    tol = 1e-6; 
    parfor j = 1:N
        f_zMp = @(x) z(j)*x - f_Mp(x);  
        [temp,~] = gmres(f_zMp,vector,[],tol,maxiter);     
        v = v+((f(z(j))/z(j))*dzdt(j))*temp;
    end
    v = (-4*K*sqrt(m*M)/(k*pi*N))*v;   
    v = imag(v);
    v = ell^(-1/p)*v;

%%%% Contour 3Nicks method 1
function [z,dzdt,K,k] = contour_points(m,M,N)
    k = (sqrt(M/m)-1)/(sqrt(M/m)+1);
    L = -log(k)/pi;
    [K,Kp] = ellipkkp(L);
    t = .5i*Kp-K+(.5:N)*2*K/N;
    [u, cn, dn] = ellipjc(t,L);
    z = sqrt(m*M)*((1/k+u)./(1/k-u));
    dzdt = cn.*dn./(1/k-u).^2;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions PCG or GMRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L_array = ichol_data(A_array,threshold)
    % Generate the Cholevsky decomposition of each matrix
    % We need this part if we want to solve A^{-k}v with a direct method
    L_array = cell(1,length(A_array));
    parfor i = 1 : length(A_array)
        L_array{i} = ichol(A_array{i},struct('type','ict','droptol',threshold,'michol','off'));
    end

function v = Mp_x_vector_pcg(A_array,L_array,vector,p)
    V = zeros(length(vector),length(A_array));
    parfor i = 1 : length(A_array)
        V(:,i) = Ap_x_vector_pcg(A_array{i},L_array{i},vector,p);
    end
    v = sum(V,2);

function v = Ap_x_vector_pcg(A,L,vector,p)
    v = vector;
    % Recall L*L' ~= A
    tol = 1e-6; maxiter = 5;
    for i = 1 : abs(p)
        [v, ~] = pcg(A,v,tol,maxiter,L,L',[]);
    end

