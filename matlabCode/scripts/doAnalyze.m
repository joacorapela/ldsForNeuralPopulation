
% PLDS toolbox example, with external input
%
% Lars Buesing, Jakob H Macke, 2014
%

clear all
close all
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animalname = 'VL63';
matfilename = 'task_remade_from_stimlaser_perf_timeSeries.mat';
binsize = '200ms/';
area = 'V1';
skipmsg = 1;
nrep = 5;

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

%%%% generate data

% trueparams = PLDSgenerateExample('xDim',xDim,'yDim',yDim,'doff',doff,'uDim',uDim);
% seqOrig    = PLDSsample(trueparams,T,Trials,'yMax',10);
% tp         = trueparams;

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

%%% check learned parameters

% figure
% plot(vec(tp.model.C*tp.model.B),vec(params.model.C*params.model.B),'xr')

%%% test model by predicting future spike trains 

% Tpred = 200;
% seqPred = PLDSsample(trueparams,Tpred,1);  % sample a test data set

% condRange = [1:100];                       % the time interval to condition on
% predRange = [101:200];                     % the time interval to predict

% predict with learned parameters
% tic; [ypred xpred xpredCov seqInf] = PLDSPredictRange(params,seqPred(1).y,condRange,predRange,'u',seqPred(1).u); toc

% figure;
% subplot(2,1,1)
% imagesc(seqPred.y(:,predRange))
% title('true data')
% subplot(2,1,2)
% imagesc(ypred)
% title('prediction')

cd(oldFolder)

clearvars -except animalname matfilename area rep skipmsg binsize
end

