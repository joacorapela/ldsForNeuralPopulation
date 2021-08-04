
% This script is called by doFit_Allanimals_trialBased which specifies
% the commented parameters (animal name, matfilename, etc.)

if ~LONO.do
    summarymatfile = ...
        dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/task_*.mat',animallist{animali},preprocessinglist{animali})));
    savedir = fullfile(resultdir,animalname,'trial_based',[num2str(binSizems),'ms','Bins']);
    savename = [area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
else
    summarymatfile = ...
        dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/LONO_*.mat',animallist{animali},preprocessinglist{animali})));
    % default: latest one 
    if length(summarymatfile) > 1
        summarymatfile = summarymatfile(end);
    end
    % add a summarymatfile name option for LONO - also save in lono folder,
    % name fold1
    savedir = fullfile(resultdir,animalname,'trial_based_LONO',[num2str(binSizems),'ms','Bins']);
    savename = ['Fold',num2str(LONO.fold),'_',area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
end


if ~isdir(savedir)
    mkdir(savedir)
end
resultsFilename = [savedir,'/',savename,'.mat'];
resultsFigname = [savedir,'/',savename,'.pdf'];


seq = buildTrialBasedSeq(summarymatfile, binSizems,binWinms,area,LONO);
config.binSizems= binSizems;
config.binWinms = binWinms;
config.area = area;
config.animalname = animalname;


codeRoot = '/mnt/data/Mitra/cache/repos/pop_spike_dyn';
oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error

uDim    = 4;
xDim    = 9;
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

params.model.inferenceHandle = @PLDSLaplaceInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

%%% save true and estimated models
if doSaveres
    save(resultsFilename, 'params', 'seq', 'varBound','config');
end
doQuickPlots(params, seq, varBound,doSavefig,resultsFigname);

cd(oldFolder)


function doQuickPlots(params, seq, varBound,doSavefig,resultsFigname)

% additional checks: check alignment of stim On per stimulus

stimOn = sum(seq(1).u(1:2,:),1); % based on trial 1

f = figure;
f.Units = 'normalized';
f.Position = [0.3988 0.0590 0.4446 0.8200];
set(f,'Color','w')
f.Name = [resultsFigname];

subplot(2,2,1);plot(varBound);
xlabel('Iteration');
ylabel('VB');

s2 = subplot(2,2,2);
plot(seq(1).posterior.xsm') % single trial
hold on;patch([find(stimOn),fliplr(find(stimOn))],[s2.YLim(1)+zeros(size(find(stimOn))),s2.YLim(2)+zeros(size(find(stimOn)))],...
    'y','FaceAlpha',0.1,'EdgeAlpha',0);
s2.Title.String = 'states variables in an example single trial';
s2.Title.FontWeight='normal';s2.Title.FontSize=9;
xlabel('timebins');

s3=subplot(2,2,3);
for state=1:9
    hold on;plot(mean(cell2mat(arrayfun(@(x) x.posterior.xsm(state,:), seq,'UniformOutput',0)'),1),'b')
    %hold on;plot(mean(cell2mat(arrayfun(@(x) nanmean(x.y,1), seq,'UniformOutput',0)'),1),'k')
end
hold on;patch([find(stimOn),fliplr(find(stimOn))],[s3.YLim(1)+zeros(size(find(stimOn))),s3.YLim(2)+zeros(size(find(stimOn)))],...
    'y','FaceAlpha',0.1,'EdgeAlpha',0);
s3.Title.String = 'state variables, averaged across all trials';
s3.Title.FontWeight='normal';s3.Title.FontSize=9;
xlabel('timebins');

s4=subplot(2,2,4);
for neuron = 1:size(seq(1).y,1)
    hold on;plot(mean(cell2mat(arrayfun(@(x) x.y(neuron,:), seq,'UniformOutput',0)'),1),'k')
end
hold on;patch([find(stimOn),fliplr(find(stimOn))],[s4.YLim(1)+zeros(size(find(stimOn))),s4.YLim(2)+zeros(size(find(stimOn)))],...
    'y','FaceAlpha',0.1,'EdgeAlpha',0);s4.Title.String = 'spiking activity, averaged across all trials';
s4.Title.FontWeight='normal';s4.Title.FontSize=9;
xlabel('timebins');

% rf0_all = (params.model.C*params.model.B);
% figure;hist(rf0_all(:,3))

if doSavefig
    saveas(f,resultsFigname,'pdf')
end
end