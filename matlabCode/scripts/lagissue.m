
clear all
%close all
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animalname = 'VL63';
matfilename = 'task_remade_from_stimlaser_perf_timeSeries.mat';
binsize = '100ms/';
area = 'V1';
skipmsg = 1;
rep = 2;% skip first one, mostly no input

% prepare 
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


% fit
% first lag 0
u4l0 = []; 
[u4l0.params,u4l0.seq,u4l0.varBound] = fitanimal(4,0,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);

u0l0 = []; 
[u0l0.params,u0l0.seq,u0l0.varBound] = fitanimal(0,0,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);

u4l1 = []; 
[u4l1.params,u4l1.seq,u4l1.varBound] = fitanimal(4,-1,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);

u0l1 = []; 
[u0l1.params,u0l1.seq,u0l1.varBound] = fitanimal(0,-1,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);

% u4l2 = []; 
% [u4l2.params,u4l2.seq,u4l2.varBound] = fitanimal(4,-2,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);
% 
% u0l2 = []; 
% [u0l2.params,u0l2.seq,u0l2.varBound] = fitanimal(0,-2,timeSeries,stateInputMemorySecs,sRate,area,trainSamples);


cd(oldFolder)
% specify the dir to save in
cd('/mnt/data/Mitra/cache/repos/figures/16072021/inputdelayissue')
save('lagissue.mat','u4l0','u0l0','u4l1','u0l1')
%save('lagissue2.mat','u4l2','u0l2')
%% plots

% load lagissue.mat here

addspikes = 0;
in = 1;
plotaligned(in,u4l0.seq,addspikes)
plotaligned(in,u4l1.seq,addspikes)
%plotaligned(in,u4l2.seq,addspikes)

figure;plot(u4l0.varBound,'b--');
hold on;plot(u0l0.varBound,'r-');
hold on;plot(u4l1.varBound,'b-');
legend({'input delay 0','no input','input delay -1'})

%%
function [params,seq,varBound] = fitanimal(un,lag,timeSeries,stateInputMemorySecs,sRate,area,trainSamples) % lag 0 or -1, un 0 or 4
    
    if strcmp(area,'V1')
        spikeCounts = [timeSeries.v1Shaft1SpikeCounts;timeSeries.v1Shaft2SpikeCounts];
    elseif strcmp(area,'LM')
        spikeCounts = [timeSeries.lmShaft1SpikeCounts;timeSeries.lmShaft2SpikeCounts];
    end

    stateInputs = buildGoNogoVisualAndLaserInputs(circshift(timeSeries.goStim,lag), ...
        circshift(timeSeries.nogoStim,lag), ...
        circshift(timeSeries.laserStim,lag), stateInputMemorySecs*sRate);

    y = spikeCounts(:,trainSamples);
    u = stateInputs(:,trainSamples);

    uDim    = un;%size(u,1);
    xDim    = 9;
    yDim    = size(y, 1);
    T       = size(y, 2);
    Trials  = 1;
    maxIter = 100;%100
    doff    = 0.0;

    seqOrig.y = y;
    seqOrig.u = u;
    seqOrig.T = T;

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
end
%%
function plotaligned(in,data,addspikes)%data=u4l1.seq
    input = data.u(in,:);
    tp = find(diff(input)>0);
    wind = [(tp'-3),(tp'-2),(tp'-1),(tp'),(tp'+1),(tp'+2),(tp'+3),(tp'+4),(tp'+5),(tp'+6)];
    figure;
    %plot(nanmean(input(wind)',2),'b');
    % this is only for vis sim
    onoff = nanmean(input(wind)',2);
    hold on;patch([find(onoff==1,1,'first'),find(onoff==1,1,'last'),find(onoff==1,1,'last'),find(onoff==1,1,'first')],...
        [-0.8 -0.8 0.8 0.8],'g','FaceAlpha',0.2,'EdgeColor','none')
    for sn=1:9
        st=data.posterior.xsm(sn,:);
        hold on;plot(nanmean(st(wind)',2)...
            ,'k');
    end

    if addspikes
        meanspike = mean(data.y,1);
        hold on;plot(nanmean(meanspike(wind)',2),'b');
    else
        figure; meanspike = mean(data.y,1);
        hold on;plot(nanmean(meanspike(wind)',2),'b.-');
        hold on;patch([find(onoff==1,1,'first'),find(onoff==1,1,'last'),find(onoff==1,1,'last'),find(onoff==1,1,'first')],...
        [-0.8 -0.8 0.8 0.8],'g','FaceAlpha',0.2,'EdgeColor','none')
 
    end
end