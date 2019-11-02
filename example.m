function example

restoredefaultpath
addpath(genpath('utils'))
addpath(genpath('subroutines'))
addpath(genpath('powerMeanLaplacian'))

%% Example where the first layer is informative and the second is not informative
numLayers                 = 2;
numNodes                  = 100;
numClusters               = 2;
sizeOfEachCluster         = numNodes/numClusters;

% Ground Truth vector
groundTruth         = [];
for j2 = 1:numClusters
    groundTruth = [groundTruth; j2*ones(sizeOfEachCluster,1)];
end
groundTruth(groundTruth == 2) = -1;

GroundTruthPerLayerCell   = {groundTruth, groundTruth};
pinVec                    = [0.9 0.1];
poutVec                   = [0.1 0.9];
sizeOfLabelSample         = 5;

s                         = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
Wcell                     = generate_multilayer_graph(numLayers, GroundTruthPerLayerCell, pinVec, poutVec);

s                         = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
idxSample                 = sample_idx_per_class(groundTruth,sizeOfLabelSample);
y                         = zeros(numNodes,1);
y(idxSample)              = groundTruth(idxSample);

% visualize adjacency matrices
figure, hold on
subplot(1,2,1), spy(Wcell{1}), title('$G^{1}$')
subplot(1,2,2), spy(Wcell{2}), title('$G^{2}$')

figure, hold on
% Clustering with p = 10
p           = 10; 
s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
C           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
test_error  = get_classification_error(C, groundTruth, idxSample)
subplot(1,5,1), stem(C), title('p=10')
1;

% Clustering with p = 1 (Arithmetic Mean)
p           = 1; 
s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
C           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
test_error  = get_classification_error(C, groundTruth, idxSample)
subplot(1,5,2), stem(C), title('p=1')
1;
1;

% Clustering with p -> 0 (log euclidean Mean)
p           = 0; 
s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
C           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
test_error  = get_classification_error(C, groundTruth, idxSample)
subplot(1,5,3), stem(C), title('p->0')
1;
1;

% Clustering with p = -1 (Harmonic Mean)
p           = -1; 
s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
C           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
test_error  = get_classification_error(C, groundTruth, idxSample)
subplot(1,5,4), stem(C), title('p=-1')
1;
1;

% Clustering with p = -10 
p           = -10; 
s           = RandStream('mcg16807','Seed',0); RandStream.setGlobalStream(s);
C           = SSL_multilayer_graphs_with_power_mean_laplacian(Wcell, p, y);
test_error  = get_classification_error(C, groundTruth, idxSample)
subplot(1,5,5), stem(C), title('p=-10')
1;
1;