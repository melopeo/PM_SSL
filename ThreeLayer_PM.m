function ThreeLayer_PM

restoredefaultpath
addpath(genpath('powerMeanLaplacian'))
addpath(genpath('subroutines'))
addpath(genpath('utils'))

dirName_Output_Data = 'ThreeLayer_PM';
if ~exist(dirName_Output_Data,'dir')
    mkdir(dirName_Output_Data)
end

% Multilayer Graph data
sizeOfEachCluster   = 100;
numClusters         = 3;
numLayers           = 3;
numNodes            = numClusters*sizeOfEachCluster;

% % % % % % % % % % % % % % % % % % % % % % % 
% GLOBAL Ground Truth vector
GroundTruth         = [];
for i = 1:numClusters
    GroundTruth = [GroundTruth; i*ones(sizeOfEachCluster,1)];
end

% Setting ground truth per layer
GroundTruthPerLayerCell     = cell(numLayers,1);
for i = 1:numLayers
    GroundTruthPerLayerCell{i} = GroundTruth == i;
end
% % % % % % % % % % % % % % % % % % % % % % % 

% Data for power means
pArray                 = [10,1,0,-1,-10];
idxNeg                 = find(pArray<=0);
lambda                 = 1;

% Setting diagonal shift depending of value of power 'p'
diagShiftArray              = zeros(size(pArray));
diagShiftArray(idxNeg)      = log10(1+abs(pArray(idxNeg)));
diagShiftArray(pArray == 0) = 1.e-6;

% Mixing parameter
diffArray         = 0:0.0025:0.1;%

pin_Layer1_Array  = (0.1+diffArray)/2; %p_in  of layer 2
pout_Layer1_Array = (0.1-diffArray)/2; %p_out of layer 2

% number of runs
numGraphRuns           = 5;
numLabelSamplesPerRun  = 5;
sizeOfLabelSampleArray = 1:0.5*sizeOfEachCluster;

for i1 = 1:length(pArray) % for per method
    p          = pArray(i1);
    diagShift  = diagShiftArray(i1);
    
    if p >= 0
        method_str = strcat( 'p_positive_', num2str(p) );
    else
        method_str = strcat( 'p_negative_', num2str(abs(p)) );
    end
    
    subdir     = strcat(dirName_Output_Data);%, '_diagShift_', num2str(diagShiftValue));
    if ~exist(subdir,'dir')
        mkdir(subdir)
    end
    
    filename_start = strcat(subdir, filesep, method_str, '_start.txt');
   
    if true%~exist(filename_start, 'file')
        save(filename_start, 'filename_start')
    
        for i2 = 1:length(diffArray) % per parameter gap
            pin  = pin_Layer1_Array(i2);
            pout = pout_Layer1_Array(i2);

            for i3 = 1:numGraphRuns % per graph run
                s        = RandStream('mcg16807','Seed',i3); RandStream.setGlobalStream(s);
                W_cell   = generate_multilayer_graph(numLayers, GroundTruthPerLayerCell, pin, pout);

                for i4 = 1:length(sizeOfLabelSampleArray) % per amount of labelled nodes
                    sizeOfLabelSample = sizeOfLabelSampleArray(i4);

                    for i5 = 1:numLabelSamplesPerRun % per run of labelled nodes
                        s                 = RandStream('mcg16807','Seed',i5); RandStream.setGlobalStream(s);
                        idxSample         = sample_idx_per_class(GroundTruth,sizeOfLabelSample);
                        y                 = zeros(numNodes,1);
                        y(idxSample)      = GroundTruth(idxSample);

                        C                 = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, y, diagShift,lambda);%, method_str, krylovOpts);

                        error_matrix(i5)  = get_classification_error(C, GroundTruth, idxSample);

                    end
                    1;
                    error_matrix_mean_labels(i4) = mean(error_matrix);    

                end
                1;
                error_matrix_mean_graph(i3,:) = error_matrix_mean_labels;
            end
            error_matrix_mean_matrix(i2,:) = mean(error_matrix_mean_graph,1);
        end
        
        filename_save = strcat(subdir, filesep, method_str, '.mat');
        save(filename_save, 'error_matrix_mean_matrix')

        createfigure(error_matrix_mean_matrix')    
        save_plots(gcf, strcat(dirName_Output_Data, filesep, method_str))

    end
    1;
end
1;

function figure1=createfigure(cdata1)

set(0,'DefaultTextInterpreter', 'latex')

% Create figure
figure1 = figure('PaperOrientation','landscape','PaperSize',[11 8.5]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
colormap jet

% Create xlabel
xL = xlabel({'$p_{\mathrm{in}} - p_{\mathrm{out}}$'});

% Create ylabel
yL = ylabel({'\% labelled nodes'});

box(axes1,'on');
axis(axes1,'tight');
axis(axes1,'square');
ax = gca;
ax.YDir = 'normal';

% Set the remaining axes properties
set(axes1, ...
    'CLim',[0 0.7], ...
    'FontSize',40, ...
    'Layer','top', ...
    'XTick',[1 21 41],...
    'XTickLabel',{'0','0.05','0.1'}, ...
    'YTick',[1 25 50], ...
    'YTickLabel',{'0', '25','50'}, ...
    'FontName', 'Helvetica');

yrule = ax.YAxis;
yrule.FontSize = 20;

xrule = ax.XAxis;
xrule.FontSize = 20;

yL.FontSize = 40;
xL.FontSize = 40;
