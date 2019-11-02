function multilayer_V1_PM_unbalanced

restoredefaultpath
addpath(genpath('powerMeanLaplacian'))
addpath(genpath('subroutines'))
addpath(genpath('utils'))

dirName_Output_Data = 'multilayer_V1_PM_unbalanced';
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

% Setting diagonal shift depending of value of power 'p'
diagShiftArray              = zeros(size(pArray));
diagShiftArray(idxNeg)      = log10(1+abs(pArray(idxNeg)));
diagShiftArray(pArray == 0) = 1.e-6;
lambda                      = 1;
loss_str                    = 'homogeneous_loss';

pin_1  = 0.09;
pout_1 = 0.01;

pin_2  = 0.05;
pout_2 = 0.05;

pin_vec  = [pin_1  pin_2];
pout_vec = [pout_1 pout_2];

% number of runs
numGraphRuns           = 10;
numLabelSamplesPerRun  = 10;
sizeOfLabelSampleArray_1 = 1:49;
sizeOfLabelSampleArray_2 = 50-sizeOfLabelSampleArray_1;

% % % % % % % % % Plot parameters % % % % % % % % % % % 
ratioLabeling = sizeOfLabelSampleArray_1./sizeOfLabelSampleArray_2;

MarkerSize       = [];
fontSize         = 30;
fontSize_legend  = 30;

xArray           = sizeOfLabelSampleArray_1;
legendLocation   = 'northoutside';
xAxisTitle_str   = 'Ratio Class Labelled Sets';%
yAxisTitle_str   = 'Classification Error';

yTickArray       = 0:0.1:0.5;
xTickArray       = [1 25 49];
xticklabels_cell = {num2str(round(100*ratioLabeling(1))/100), '1',num2str(round(100*ratioLabeling(end))/100)};

title_str        = '';
legend_boolean   = true;

[legendCell, colorCell, markerCell, LineStyleCell, LineWidthCell] =get_plot_parameters_SPM;
modelsToignore                 = [2 6];
legendCell(modelsToignore)     = [];
colorCell(modelsToignore)      = [];
markerCell(modelsToignore)     = [];
LineStyleCell(modelsToignore)  = [];
LineWidthCell(modelsToignore)  = [];


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

filename_start = strcat(dirName_Output_Data, filesep, 'start.txt');
if true%~exist(filename_start, 'file')
	save(filename_start, 'filename_start')

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
	   
	%    if ~exist(filename_start, 'file')
	%        save(filename_start, 'filename_start')
		
		for i2 = 1:length(sizeOfLabelSampleArray_1)
		    sizeOfLabelSample(1) = sizeOfLabelSampleArray_1(i2);

		        i4 = 1;
		        sizeOfLabelSample(2) = sizeOfLabelSampleArray_2(i2);

		        for i3 = 1:numGraphRuns
		            s        = RandStream('mcg16807','Seed',i3); RandStream.setGlobalStream(s);
		            W_cell   = generate_multilayer_graph(numLayers, GroundTruthPerLayerCell, pin_vec, pout_vec);

		            for i5 = 1:numLabelSamplesPerRun
		                s                 = RandStream('mcg16807','Seed',i5); RandStream.setGlobalStream(s);
		                idxSample         = sample_idx_per_class(GroundTruth,sizeOfLabelSample);
		                y                 = zeros(numNodes,1);
		                y(idxSample)      = GroundTruth(idxSample);

		                C                 = SSL_multilayer_graphs_with_power_mean_laplacian(W_cell, p, y, diagShift, lambda, loss_str);

		                error_matrix_mean_labels(i3,i5) = get_classification_error(C, GroundTruth, idxSample);
		            end

		        1;
		        error_matrix_mean_graph(i2,i4) = mean(error_matrix_mean_labels(:));
		    end
		end
		error_matrix_mean_graph_matrix(:,i1) = error_matrix_mean_graph;
		1;
	%    end
	    1;
	end
	1;
	% Plot
	fig_handle      = plot_performance(error_matrix_mean_graph_matrix', xArray, legendCell, colorCell, LineStyleCell, markerCell, MarkerSize, LineWidthCell,legendLocation,xAxisTitle_str,yAxisTitle_str, title_str, fontSize,fontSize_legend,legend_boolean,xTickArray,yTickArray,xticklabels_cell);
	filename_prefix = strcat(dirName_Output_Data, filesep, 'output');
	save_plots(fig_handle, filename_prefix)
end

function fig_handle = plot_performance(mean_Matrix, xArray, legendCell, colorCell, LineStyleCell, markerCell, MarkerSize, LineWidthCell,legendLocation,xAxisTitle_str,yAxisTitle_str, title_str, fontSize,fontSize_legend,legendBoolean,xTickArray,yTickArray,xticklabels_cell)

    fig_handle = figure; hold on

    for j = 1:size(mean_Matrix,1)
        meanVec = mean_Matrix(j,:);

        plot(xArray,meanVec, ...
            'Color',colorCell{j}, ...
            'Marker', markerCell{j}, ...
            'MarkerFaceColor', colorCell{j}, ...
            'MarkerEdgeColor',colorCell{j}, ...
            'LineWidth', LineWidthCell{j}, ...
            'LineStyle', LineStyleCell{j}, ...
            'MarkerSize',10);
    end
    set(gca,'XTick', xArray(xTickArray))
    xticklabels(xticklabels_cell)
    if legendBoolean
        legend(legendCell,'Location',legendLocation, 'Interpreter','latex','FontSize',fontSize_legend, 'Orientation', 'horizontal', 'Location', 'northoutside')
    end
    axis square tight
    box on
    daspect

    ax = gca;
    xlabel(xAxisTitle_str, 'interpreter', 'latex', 'FontSize', fontSize,'fontweight','bold')
    ylabel(yAxisTitle_str, 'interpreter', 'latex', 'FontSize', fontSize,'fontweight','bold')
    title(title_str, 'interpreter', 'latex', 'FontSize', fontSize)
    set(gca,'fontweight','bold', 'FontSize', fontSize);
