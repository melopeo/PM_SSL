function multilayer_V1_PM

restoredefaultpath
addpath(genpath('powerMeanLaplacian'))
addpath(genpath('subroutines'))
addpath(genpath('utils'))

dirName_Output_Data = 'multilayer_V1_PM';
if ~exist(dirName_Output_Data,'dir')
    mkdir(dirName_Output_Data)
end

% Multilayer Graph data
sizeOfEachCluster   = 100;
numClusters         = 2;
numLayers           = 2;
numNodes            = numClusters*sizeOfEachCluster;

% Ground Truth vector
GroundTruth         = [];
for j2 = 1:numClusters
    GroundTruth = [GroundTruth; j2*ones(sizeOfEachCluster,1)];
end
GroundTruth(GroundTruth == 2) = -1;

% Setting ground truth per layer
GroundTruthPerLayerCell     = cell(numLayers,1);
for j2 = 1:numLayers
    GroundTruthPerLayerCell{j2} = GroundTruth;
end

% Data for power means
pArray                 = [10,1,0,-1,-10];
idxNeg                 = find(pArray<=0);
lambda                 = 1;

% Setting diagonal shift depending of value of power 'p'
diagShiftArray              = zeros(size(pArray));
diagShiftArray(idxNeg)      = log10(1+abs(pArray(idxNeg)));
diagShiftArray(pArray == 0) = 1.e-6;

% Mixing parameter
diffArray         = -0.1:0.005:0.1;

pin_Layer1_Array  = (0.1+diffArray)/2; %p_in  of layer 2
pout_Layer1_Array = (0.1-diffArray)/2; %p_out of layer 2

pin_Layer2_Array  = 0.09*ones(size(diffArray)); %p_in  of layer 1
pout_Layer2_Array = 0.01*ones(size(diffArray)); %p_out of layer 1

% number of runs
numGraphRuns           = 10;
numLabelSamplesPerRun  = 10;
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
            pinVec  = [pin_Layer1_Array(i2)  pin_Layer2_Array(i2)];
            poutVec = [pout_Layer1_Array(i2) pout_Layer2_Array(i2)];

            for i3 = 1:numGraphRuns % per grapn run
                s        = RandStream('mcg16807','Seed',i3); RandStream.setGlobalStream(s);
                W_cell   = generate_multilayer_graph(numLayers, GroundTruthPerLayerCell, pinVec, poutVec);

                for i4 = 1:length(sizeOfLabelSampleArray) % per amount of labelled nodes
                    sizeOfLabelSample = sizeOfLabelSampleArray(i4);

                    for i5 = 1:numLabelSamplesPerRun % per run of labelled nodes
                        s                 = RandStream('mcg16807','Seed',i5); RandStream.setGlobalStream(s);
                        idxSample         = sample_idx_per_class(GroundTruth,sizeOfLabelSample);
                        y                 = zeros(numNodes,1);
                        y(idxSample)      = GroundTruth(idxSample);

                        % classify
                        C                 = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, y, diagShift, lambda);%, method_str, krylovOpts);

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
        
        figure_handle=createfigure(error_matrix_mean_matrix');
        save_plots(figure_handle, strcat(dirName_Output_Data, filesep, method_str))
        1;
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
xL = xlabel({'$G^{(2)}:$\,\, $p_{\mathrm{in}}^{(2)} - p_{\mathrm{out}}^{(2)}$'});

% Create ylabel
yL = ylabel({'\% labelled nodes'});

box(axes1,'on');
axis(axes1,'tight');
axis(axes1,'square');
ax = gca;
ax.YDir = 'normal';

% Set the remaining axes properties
set(axes1, ...
    'CLim',[0 0.5], ...
    'FontSize',40, ...
    'Layer','top', ...
    'XTick',[1 21 41],...
    'XTickLabel',{'-0.1','0','0.1'}, ...
    'YTick',[1 25 50], ...
    'YTickLabel',{'0', '25','50'}, ...
    'FontName', 'Helvetica');

yrule = ax.YAxis;
yrule.FontSize = 20;

xrule = ax.XAxis;
xrule.FontSize = 20;

yL.FontSize = 40;
xL.FontSize = 40;

% title
title_str = '\fontsize{15}{0} $ G^{(1)}: p_{\mathrm{in}}^{(1)}=0.09,\,\,  p_{\mathrm{out}}^{(1)}=0.01 $';
title(title_str, 'interpreter', 'latex', 'FontSize', 22)

% add lines
vline(21,'k-')
