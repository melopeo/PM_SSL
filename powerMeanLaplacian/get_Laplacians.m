function matrix_cell = get_Laplacians(adjacency_cell, diagShift)
% matrix_cell = get_Laplacians(adjacency_cell, str_cell, p)
% This function gets the Laplacian or Signed Laplacian for each adjacency
% matrix given in adjacency_cell, and adds a shift according to the
% parameters p
% INPUT : adjacency_cell :  cell array of adjacency matrices
%       : str_cell       :  cell array of strings: entries stand for
%                                 : 'L': for Laplacian
%       : p              : integer

if nargin < 2
    diagShift = 0;
end

numMatrices = length(adjacency_cell);          % number of layers
n           = size(adjacency_cell{1},1);       % number of nodes
matrix_cell = cell(numMatrices,1);

for i = 1:numMatrices
    
    W              = adjacency_cell{i};         % adjacency matrix
    L              = build_laplacian_matrix(W); % Laplacian matrix
    matrix_cell{i} = L;

end

% apply shift
matrix_cell = cellfun(@(x) x + diagShift*speye(n), matrix_cell, 'UniformOutput', false);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function L = build_laplacian_matrix(W)
% L = build_Laplacian_matrix(W)
% This function builds a sparse Laplacian Matrix
% input:  W (matrix)    : adjacency matrix
% output: L (matrix)    : sparse Laplacian matrix
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

n            = size(W,1);              % size of graph
d            = sum(W,2);               % degree vector
dInv         = 1./d;                   % d^{-1}                   
dInv(d == 0) = 0;                      % take care of isolated nodes

dInv         = dInv.^(0.5);            % d^{-1/2}
DInv         = spdiags(dInv, 0, n, n); % D^{-1/2}

% Symmetric Laplacian
L = spdiags(d ~= 0, 0, n, n) - DInv*W*DInv; 

% enforce symmetry
L = (L  + L')/2;
