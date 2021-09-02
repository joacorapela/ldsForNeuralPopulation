function [params ,seq ,varBound ,EStepTimes ,MStepTimes] = dofitWithNstates(xDim,seq)

uDim    = size(seq(1).u,1);
% xDim    = 9;
yDim    = size(seq(1).y,1);
T       = size(seq(1).T);
Trials  = length(seq);
maxIter = 100;
doff    = 0.0;

fprintf('Max spike count:    %i \n', max(vec([seq.y])))
fprintf('Mean spike count:   %d \n', mean(vec([seq.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([seq.y])>0.5))

%%% fit model

params = [];
% important: set flag to use external input
if uDim>0;params.model.notes.useB = true;end

params = PLDSInitialize(seq, xDim, 'NucNormMin', params);
% fprintf('Initial subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

params.model.inferenceHandle = @PLDSVariationalInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))