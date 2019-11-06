function realworld_experiments

restoredefaultpath
addpath(genpath('powerMeanLaplacian'))
addpath(genpath('subroutines'))


% datasets location and info
dir_data2               = 'realworld_datasets';
dataname_cell           = {'3sources','BBC4view_685','BBCSport2view_544','WikipediaArticles', 'UCI_mfeat', 'citeseer', 'cora', 'webKB_texas_2'};
dataname_cell_for_print = {'3sources','BBC','BBCS','WikipediaArticles', 'UCI', 'Citeseer', 'Cora', 'WebKB'};

% general settings
knn                    = 10;
numSampleRuns          = 10;
sizeOfLabelSampleArray = [0.01 0.05:0.05:0.25];

% Data for power means
pArray                 = [1,-1,-10];
idxNeg                 = find(pArray<=0);
lambda_array           = [0.1 10 10];

% Setting diagonal shift depending of value of power 'p'
diagShiftArray              = zeros(size(pArray));
diagShiftArray(idxNeg)      = log10(1+abs(pArray(idxNeg)));
diagShiftArray(pArray == 0) = 1.e-6;

formatSpec = 'Dataset: %s - Labeled Nodes: %3.0f %% - Power(p): %d - Average Classification error: %3.1f %% \n';

for i1 = 1:length(dataname_cell) % per dataset
    
    dataname                = dataname_cell{i1};
    dataname_for_print      = dataname_cell_for_print{i1};  
        
    dataname_file = strcat(dir_data2, filesep, dataname, filesep, 'knn_', num2str(knn), '.mat');

    dataset       = load(dataname_file);
    W_cell        = dataset.W_cell;
    labels        = dataset.labels;
    numNodes      = size(W_cell{1},1);
 
       for i2 = 1:length(pArray) % per power

           p            = pArray(i2);
           diagShift    = diagShiftArray(i2);
           lambda       = lambda_array(i2);
                
            for i3 = 1:length(sizeOfLabelSampleArray) % per training data

                sizeOfLabelSample = sizeOfLabelSampleArray(i3);
                error_C           = inf(numSampleRuns,1);
                
                for i4 = 1:numSampleRuns % per run of labeled nodes

                    s = RandStream('mcg16807','Seed',i4); RandStream.setGlobalStream(s);
                    idxSample         = sample_idx_per_class(labels, sizeOfLabelSample, 'percentage');
                    y                 = zeros(numNodes,1);
                    y(idxSample)      = labels(idxSample);
                    1;

                    s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
                    C           = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, y, diagShift, lambda);
                    error_C(i4) = get_classification_error(C, labels(:), idxSample);

                end

                mean_error_C = mean(error_C);

                fprintf(formatSpec, dataname_for_print, 100*sizeOfLabelSample, p, 100*mean_error_C)
                1;
            end 
            1;
            fprintf('\n')
       end
end
