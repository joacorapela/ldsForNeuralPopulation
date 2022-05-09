
% This script is called by doFit_Allanimals_trialBased which specifies
% the commented parameters (animal name, matfilename, etc.)

summarymatfile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/task_*.mat',animallist{animali},preprocessinglist{animali})));
% only for correct trials
 perffile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/rev_perf*.mat',animallist{animali},preprocessinglist{animali})));


if ~CntrlOnly
    savedir = fullfile(resultdir,animalname,['Joint_trial_based_splitContext'],[num2str(binSizems),'ms','Bins']);
else
    savedir = fullfile(resultdir,animalname,['Joint_trial_based_splitContext_CntrlOnly'],[num2str(binSizems),'ms','Bins']);
end
savename = ['Joint_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];


if ~isdir(savedir)
    mkdir(savedir)
end

splitDelays = 0;
% if Xval.do ==1 seq will be an Xval struct

[seq,N_V1,N_LM] = buildTrialBasedSeq_JointSepCon(summarymatfile, binSizems,binWinms,splitDelays,perffile,trialType,Xval,CntrlOnly);

savename = [savename,'_RND',num2str(RandSeed)];
savename = [savename,'_',trialType];

config.binSizems= binSizems;
config.binWinms = binWinms;
% config.area = area;
config.animalname = animalname;
config.summarymatfile = summarymatfile;
config.splitDelays = splitDelays;
config.baselineU = baselineU;
config.trialType = trialType;
config.RandSeed = RandSeed;
config.CntrlOnly = CntrlOnly;

codeRoot = '/mnt/data/Mitra/cache/repos/pop_spike_dyn';
oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error

%if numel(nStates) == 1 % evaluatig single model

if Xval.do
    
    resultsFilename = [savedir,'/',savename,'_',num2str(Xval.nFolds),'FoldXval','.mat'];

    Xval = seq;
    Xval.params = cell(1,Xval.nFolds);
    Xval.varBound = cell(1,Xval.nFolds);
     for FoldN = 1:Xval.nFolds
         resultsFigname = [savedir,'/',savename,'_XvalFold',num2str(FoldN),'.pdf'];
         
         [Xval.params{FoldN} ,Xval.train_seq{FoldN},Xval.varBound{FoldN} ,EStepTimes ,MStepTimes] = ...
             dofitWithNstates_JointSepCon(nStates,Xval.train_seq{FoldN},Inference_handle,config,N_V1,N_LM);
       
        try 
            doQuickPlots(Xval.params{FoldN}, Xval.train_seq{FoldN}, Xval.varBound{FoldN},doSavefig,resultsFigname,CntrlOnly);
        catch
        end
         
     end
     if doSaveres
         save(resultsFilename, 'Xval','config');
     end
    
else
    resultsFilename = [savedir,'/',savename,'.mat'];
    resultsFigname = [savedir,'/',savename,'.pdf'];


    [params ,seq ,varBound ,EStepTimes ,MStepTimes] = dofitWithNstates_JointSepCon(nStates,seq,Inference_handle,config,N_V1,N_LM);
    %%% save true and estimated models
    if doSaveres
        save(resultsFilename, 'params', 'seq', 'varBound','config');
    end
    doQuickPlots(params, seq, varBound,doSavefig,resultsFigname,CntrlOnly);
end
% else % model selection
%     Allmodels = cell(1,numel(nStates));
%     count = 1;
%     for nst = nStates
%         clear resTr
%         FittedFold = LONO.fold;
%         try % it might error with some nsts
%             [resTr.params ,resTr.seq ,resTr.varBound ,resTr.EStepTimes ,resTr.MStepTimes] = dofitWithNstates(nst,seq,Inference_handle,config);
%             resTr.LONO = LONO;
%             resTr.config = config;
%             doLONO;
%             %doLONO_and_ValLogLik;
%             resTr.lono_trial_ll = trial_ll;
%             resTr.lono_Avtrial_ll = mean(mean(trial_ll,2));
%             resTr.test_trial_ll = test_trial_ll;
%             resTr.test_Avtrial_ll = mean(mean(test_trial_ll,2));
%         catch
%         end
%         resTr.nStates = nst;
%         % save trial_ll along with ns and some state variables
%         % also add option for repeating folds
%         Allmodels{count} = resTr;
%         count = count + 1;
%     end
%     try 
%         save(resultsFilename,'Allmodels')
%     catch
%         save(resultsFilename,'Allmodels','-v7.3')
%     end
% end

cd(oldFolder)
function doQuickPlots(params, seq, varBound,doSavefig,resultsFigname,CntrlOnly)

% additional checks: check alignment of stim On per stimulus
if ~CntrlOnly
    stimOn = sum(seq(1).u(1:2,:),1); % based on trial 1
else
    stimOn = sum(seq(1).u(1,:),1); % based on trial 1
end

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