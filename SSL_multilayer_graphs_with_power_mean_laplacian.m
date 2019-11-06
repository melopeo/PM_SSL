function f = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, labels_vec, diagShift, lambda, loss_str)
% C = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, labels_vec, diagShift, lambda, post_processing_str)
% INPUT: W_cell (cell)              : cell containing adjacency matrices
%        p (scalar)                 : p-th power of generalized matrix mean
%        labels_vec (array)         : array with node labels (zero entries for unlabelled nodes)
%        diagShift (scalar)         : diagonal shift of Laplacians for the case where p<0
%                                     deafult:
%                                     - for p <= 0 : log10(1+abs(p))+1.e-6
%                                     - for p > 0  : zero
%        lambda (scalar)            : regularizer parameter (default value lambda = 1)
%                                       eigenvectors of power mean Laplacian Lp
% OUTPUT: f (array)                 : Classification output                            

% Reference:

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % Process input parameters % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% if diagonal shift is not given, set it depending on value of p
if nargin < 4 
    if p > 0
        diagShift = 0;
    else
        diagShift = log10(1+abs(p))+1.e-6;
    end
end

if nargin < 5
    lambda = 1;
end

if nargin < 6
    loss_str = 'homogeneous_loss';
end

classes_idx = unique(labels_vec(labels_vec ~= 0));
num_classes = length(classes_idx);

n           = size(W_cell{1},1);                       % nmber of nodes
matrix_cell = get_Laplacians(W_cell, diagShift);       % compute Laplacians
M           = get_matrix_power_mean(matrix_cell, p);   % get matrix power mean
   
% classification: one versus all
F = zeros(n,num_classes);

if strcmp('homogeneous_loss', loss_str)
    for i = 1:num_classes % classification per class
        label_vec_class_i = (labels_vec == classes_idx(i));
        F(:,i) = (eye(n)+lambda*M)\label_vec_class_i;       % (I+lambda*M)*f = y
    end
    [~, f_aux] = max(F,[],2);
    
elseif strcmp('heterogeneous_loss', loss_str)
    for i = 1:num_classes % classification per class
        label_vec_class_i = (labels_vec == classes_idx(i))+0;
        label_vec_class_i(label_vec_class_i==1) = n/sum(label_vec_class_i); % weighted class vector
        F(:,i) = (eye(n)+lambda*M)\label_vec_class_i;       % (I+lambda*M)*f = y
    end
    [~, f_aux] = max(F,[],2);
    
elseif strcmp('CMN', loss_str)
    for i = 1:num_classes
        label_vec_class_i = (labels_vec == classes_idx(i))+0;
        fl(:,i)           = label_vec_class_i;
        F(:,i) = (eye(n)+lambda*M)\label_vec_class_i;       % (I+lambda*M)*f = y
    end
    
    % class mass normalization
    q = sum(fl)+1; % the unnormalized class proportion estimate from labeled data, with Laplace smoothing
    fu = F(labels_vec==0,:);
    l  = sum(labels_vec~=0); % the number of labeled points
    fu_CMN = fu .* repmat(q./sum(fu), n-l, 1);
    
    [~, f_aux] = max(fu_CMN,[],2);

end
    
% % % % % % Formatting output: relabel
f = zeros(n,1);
for i = 1:num_classes
    loc    = f_aux == i;
    f(loc) = classes_idx(i);
end

