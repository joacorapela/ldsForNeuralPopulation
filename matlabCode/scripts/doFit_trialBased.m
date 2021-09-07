
% This script is called by doFit_Allanimals_trialBased which specifies
% the commented parameters (animal name, matfilename, etc.)

summarymatfile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/task_*.mat',animallist{animali},preprocessinglist{animali})));

if ~LONO.do
    savedir = fullfile(resultdir,animalname,['trial_based_split_',num2str(splitDelays)],[num2str(binSizems),'ms','Bins']);
    savename = [area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
else
    LONO.file = ...
        dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/LONO_*.mat',animallist{animali},preprocessinglist{animali})));
    % default: latest one
    if length(LONO.file) > 1
        LONO.file = LONO.file(end);
    end
    % add a summarymatfile name option for LONO - also save in lono folder,
    % name fold1
    savedir = fullfile(resultdir,animalname,['trial_based_MS_split_',num2str(splitDelays)],[num2str(binSizems),'ms','Bins']);
    savename = ['Fold',num2str(LONO.fold),'_',area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
    if numel(nStates) > 1 % doing model selection
        savename = ['ModelSelection-',savename];
    end
end


if ~isdir(savedir)
    mkdir(savedir)
end
resultsFilename = [savedir,'/',savename,'.mat'];
resultsFigname = [savedir,'/',savename,'.pdf'];


seq = buildTrialBasedSeq(summarymatfile, binSizems,binWinms,area,LONO,splitDelays);
config.binSizems= binSizems;
config.binWinms = binWinms;
config.area = area;
config.animalname = animalname;
config.summarymatfile = summarymatfile;
config.splitDelays = splitDelays;

codeRoot = '/mnt/data/Mitra/cache/repos/pop_spike_dyn';
oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error

if numel(nStates) == 1 % evaluatig single model
    [params ,seq ,varBound ,EStepTimes ,MStepTimes] = dofitWithNstates(nStates,seq,Inference_handle);
    %%% save true and estimated models
    if doSaveres
        save(resultsFilename, 'params', 'seq', 'varBound','config','LONO');
    end
    doQuickPlots(params, seq, varBound,doSavefig,resultsFigname);
    
else % model selection
    Allmodels = cell(1,numel(nStates));
    count = 1;
    for nst = nStates
        clear resTr
        FittedFold = LONO.fold;
        try % it might error with some nsts
            [resTr.params ,resTr.seq ,resTr.varBound ,resTr.EStepTimes ,resTr.MStepTimes] = dofitWithNstates(nst,seq,Inference_handle);
            resTr.LONO = LONO;
            resTr.config = config;
            doLONO_and_ValLogLik;
            resTr.lono_trial_ll = trial_ll;
            resTr.lono_Avtrial_ll = mean(mean(trial_ll,2));
            resTr.test_trial_ll = test_trial_ll;
            resTr.test_Avtrial_ll = mean(mean(test_trial_ll,2));
        catch
        end
        resTr.nStates = nst;
        % save trial_ll along with ns and some state variables
        % also add option for repeating folds
        Allmodels{count} = resTr;
        count = count + 1;
    end
    try 
        save(resultsFilename,'Allmodels')
    catch
        save(resultsFilename,'Allmodels','-v7.3')
    end
end

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