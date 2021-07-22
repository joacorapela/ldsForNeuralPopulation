
% This script is usually called by doAnalyzeAllAnimals.m which specifies
% the commented parameters (animal name, matfilename, etc.)

clear all
%close all
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animalname = 'VL63';
matfilename = 'task_remade_from_stimlaser_perf_timeSeries.mat';
binsize = '100ms/';
area = 'V1';
skipmsg = 1;
nrep = 2;% skip first one, mostly no input

for rep = 2:nrep

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

stateInputs = buildGoNogoVisualAndLaserInputs(circshift(timeSeries.goStim,-1), ...
    circshift(timeSeries.nogoStim,-1), ...
    circshift(timeSeries.laserStim,-1), stateInputMemorySecs*sRate);

y = spikeCounts(:,trainSamples);
u = stateInputs(:,trainSamples);

uDim    = 0;%size(u,1);
xDim    = 9;
yDim    = size(y, 1);
T       = size(y, 2);
Trials  = 1;
maxIter = 25;%100
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
% 
% y = seqOrig.y;
% u = seqOrig.u;
% save(resultsFilename, 'y', 'u', 'params', 'seq', 'varBound');


cd(oldFolder)

%clearvars -except animalname matfilename area rep nrep skipmsg binsize matfilenamelist animallist animali
end

figure;plot(varBound)
%%
figure;plot(u');
hold on;plot(mean(y,1))

figure;plot(u');
for i = 1:9
    hold on;plot(seq.posterior.xsm(i,:))
end

%% do a model based psth
in = 1;
input = u(in,:);
tp = find(diff(input)>0);
wind = [(tp'-3),(tp'-2),(tp'-1),(tp'),(tp'+1),(tp'+2),(tp'+3),(tp'+4),(tp'+5),(tp'+6)];
figure;plot(nanmean(input(wind)',2));
for sn=1:9
st=seq.posterior.xsm(sn,:);
hold on;plot(nanmean(st(wind)',2));
end

meanspike = mean(y,1);
hold on;plot(nanmean(meanspike(wind)',2),'k');

%%%% binned data is misaligned!
%% look at binning
figure;plot(nanmean(timeSeries.v1Shaft2SpikeCounts(:,5000:10000),1))
hold on;plot(circshift(timeSeries.goStim(5000:10000),-1))
