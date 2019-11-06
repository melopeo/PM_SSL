function idxSample = sample_idx_per_class(GroundTruth,sizeOfLabelSample, style_str)

% possible vlaues for style_str: 'original', 'percentage'
% if original: then sample 'sizeOfLabelSample' elements per class
% percentage: then sample a 'sizeOfLabelSample' percentage of elements per
% class

% process inputs
if nargin < 3
    style_str = 'original';
end

% if isscalar(sizeOfLabelSample)
%     classes_idx       = unique(GroundTruth);
%     classes_num       = length(classes_idx);
%     sizeOfLabelSample = sizeOfLabelSample*ones(classes_num,1);
% end

% start
classes_idx = unique(GroundTruth);
classes_num = length(classes_idx);

if strcmp(style_str, 'original')
    
    if isscalar(sizeOfLabelSample)
        classes_idx       = unique(GroundTruth);
        classes_num       = length(classes_idx);
        sizeOfLabelSample = sizeOfLabelSample*ones(classes_num,1);
    end
  
    for i = 1:classes_num
        label_i           = classes_idx(i);
        classes_idx_i     = find( GroundTruth == label_i );
        idxSample_cell{i} = randsample(classes_idx_i, sizeOfLabelSample(i));
    end

%     idxSample = cell2mat(idxSample_cell);
%     idxSample = sort(idxSample, 'ascend');
%     idxSample = idxSample(:);
    
else
    for i = 1:classes_num
        label_i             = classes_idx(i);
        classes_idx_i       = find( GroundTruth == label_i );
        class_size_i        = length(classes_idx_i);
        sizeOfLabelSample_i = floor(sizeOfLabelSample*class_size_i);
        % we want sizeOfLabelSample_i to be smaller than class size and larger than 0
        sizeOfLabelSample_i = min(sizeOfLabelSample_i, class_size_i);
        sizeOfLabelSample_i = max(sizeOfLabelSample_i, 1);
        idxSample_cell{i}   = randsample(classes_idx_i, sizeOfLabelSample_i);
    end
    
end

try
    idxSample = cell2mat(idxSample_cell);
catch
    idxSample = cell2mat(idxSample_cell');
end

idxSample = sort(idxSample, 'ascend');
idxSample = idxSample(:);