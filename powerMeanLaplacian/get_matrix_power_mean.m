function M = get_matrix_power_mean(matrix_cell, p)
% [V,D] = get_power_mean_matrix(matrix_cell, p, k, method_str, krylovOpts)
% This function computes the eigenvectores corresponding to the k smallest
% eigenvalues of the generalized matrix mean (0.5*(A^p + B^p))^(1/p)
% INPUT:  matrix_cell (cell)  : cell array (k,1) in each entry it contains a
%                             symmetric positive semidefinite matrix A_i
%         p (scalar)          : p-th power of generalized matrix mean
% OUTPUT: M                   : matrix power mean

M = avg_matrix_power(matrix_cell, p);

if p ~= 0
    M = power_of_a_matrix(full(M),1/p);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function M = avg_matrix_power(matrix_cell, p)
% M = average_matrix_power(matrix_cell, p, tol)
% This function approximates the average or power matrices, i.e.
%            x = (1/k)*[(A_1)^p + ... +(A_k)^p]
% INPUT: matrix_cell : cell array (k,1) in each entry it contains a
%                      symmetric positive semidefinite matrix A_i
%        p           : scalar [power parameter]
% 

if p == 0
    % compute log euclidean mean
    M = get_matrix_exp_log_mean(matrix_cell);
    
else

    % compute average matrix power
    numberOfMatrices = length(matrix_cell);    % number of matrices in cell
    n                = size(matrix_cell{1},1); % size of matrices
    M                = zeros(n,n);
    for i = 1:numberOfMatrices
        A = matrix_cell{i};
        M = M + power_of_a_matrix(full(A),p);
    end
    M = 0.5*(M+M');                            % enforcing symmetry
    M = M/numberOfMatrices;
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function M = get_matrix_exp_log_mean(matrix_cell)
% M = get_matrix_exp_log_mean(matrix_cell)
% This function approximates the log euclidean mean, i.e.
%            x = exp( (1/k)*[log(A_1) + ... +log(A_k)] )
% INPUT: matrix_cell : cell array (k,1) in each entry it contains a
%                      symmetric positive definite matrix A_i
% 

numberOfMatrices = length(matrix_cell);    % number of matrices in cell
n                = size(matrix_cell{1},1); % size of matrices
M                = zeros(n,n);
for i = 1:numberOfMatrices
    A     = matrix_cell{i};
    logmA = logm(full(A));
    logmA = 0.5*(logmA+logmA');            % enforcing symmetry
    M     = M + logmA;
end
M = M/numberOfMatrices;
M = expm(M);
M = 0.5*(M+M');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function FH = power_of_a_matrix(H,p)
% FH = power_of_a_matrix(H,p)
% This function calculates the power of a matrix, H^p, based on the spectral
% decomposition. 
% INPUT:  H   : matrix ( n x n )
%         p   : scalar
% OUTPUT: FH  : matrix ( n x n )

[V,D] = eig(H, 'vector');
Dp    = D.^p;
FH    = V*diag(Dp)*V';
FH    = 0.5*(FH + FH');

