
% This script is usually called by doAnalyzeAllAnimals.m which specifies
% the commented parameters (animal name, matfilename, etc.)


for rep = 1:nrep

stateInputMemorySecs = 0;
trainDurSecs = 180; %180,300,600
analysisStartTimeSecs = (rep)*trainDurSecs;  % or (rep-1), from 0
if  ~isdir(['../../results/',animalname,'/',binsize,num2str(trainDurSecs)])
    mkdir(['../../results/',animalname,'/',binsize,num2str(trainDurSecs)])
end

codeRoot = '/mnt/data/Mitra/cache/repos/pop_spike_dyn';
timeSeriesFilename = ['../../results/',animalname,'/',binsize,matfilename];
resultsFilename = ['../../results/',animalname,'/',binsize,num2str(trainDurSecs),'/',area,'_PLDSresults_',num2str(analysisStartTimeSecs),'.mat'];

timeSeries = load(timeSeriesFilename);
timeSeries = timeSeries.timeSeries;
sRate = timeSeries.sRate;
trainSamples = (analysisStartTimeSecs*sRate)+(1:(trainDurSecs*sRate));

oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error

if strcmp(area,'V1')
    spikeCounts = [timeSeries.v1Shaft1SpikeCounts;timeSeries.v1Shaft2SpikeCounts];
elseif strcmp(area,'LM')
    spikeCounts = [timeSeries.lmShaft1SpikeCounts;timeSeries.lmShaft2SpikeCounts];%timeSeries.v1Shaft1SpikeCounts;
end

stateInputs = buildGoNogoVisualAndLaserInputs(timeSeries.goStim, timeSeries.nogoStim, timeSeries.laserStim, stateInputMemorySecs*sRate);

y = spikeCounts(:,trainSamples);
u = stateInputs(:,trainSamples);

uDim    = size(u,1);
xDim    = 9;
yDim    = size(y, 1);
T       = size(y, 2);
Trials  = 1;
maxIter = 100;
doff    = 0.0;

seqOrig.y = y;
seqOrig.u = u;
seqOrig.T = T;

fprintf('Max spike count:    %i \n', max(vec([seqOrig.y])))
fprintf('Mean spike count:   %d \n', mean(vec([seqOrig.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([seqOrig.y])>0.5))

%%% fit model

seq    = seqOrig;
params = [];
% important: set flag to use external input
if uDim>0;params.model.notes.useB = true;end

params = PLDSInitialize(seq, xDim, 'NucNormMin', params);
% fprintf('Initial subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

params.model.inferenceHandle = @PLDSLaplaceInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

%%% save true and estimated models

y = seqOrig.y;
u = seqOrig.u;
save(resultsFilename, 'y', 'u', 'params', 'seq', 'varBound');


cd(oldFolder)
clearvars -except animalname matfilename area rep nrep skipmsg binsize matfilenamelist animallist animali
end

