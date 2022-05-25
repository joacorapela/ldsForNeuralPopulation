animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

rootdir = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/LONO_LOAO';
%% 0 - plotting model prediction and compare with data trial-averages, for both LONO and LOAO

%%% params
animali = 1;
leaveout = 'neuron'; % neuron or area: all V1
CntrlOnly = 0; % 0 or not
tSilencingLag = 2;
onsetinbin = 62; % for lag2
dur = 9; % approximately
includeslc = 0;

%%% load the lono/loao file here
filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
res = load(fullfile(rootdir,filename));
Ypred = res.Ypred;
TrYcntrl = res.TrYcntrl;
% this is dependent on tSilencingLag
Ypred_slc = res.Lag{tSilencingLag}.Ypred_slc;
TrY = res.Lag{tSilencingLag}.TrY;

figure;
NeuronList = [5 6 11 14 15 18 19 20 23 25 26];%[5 11 12 15 18 19 20 23 24 25];% 34 for animal 6, 22 qnd 24 for animal 3
for NeuronNum = 1:10
    % trial averaged
    subplot(10,1,NeuronNum);
    plot(nanmean(Ypred(:,:,NeuronList(NeuronNum)),1),'color',[.7 .7 .7])
    hold on;plot(nanmean(TrYcntrl(:,:,NeuronList(NeuronNum)),1),'k')
    if includeslc
        hold on;plot(nanmean(Ypred_slc(:,:,NeuronList(NeuronNum)),1),'color','c')
        hold on;plot(nanmean(TrY(:,:,NeuronList(NeuronNum)),1),'b')
    end
    % this is only relevant for t2, change
    yl = get(gca,'YLim');
    if includeslc
        hold on;patch([onsetinbin,onsetinbin+dur,onsetinbin+dur,onsetinbin],[yl(1) yl(1) yl(2) yl(2)],[0 .8 0.8],'FaceAlpha',.2,'EdgeAlpha',0);
    end
    hold on;patch([116/2,3*116/4,3*116/4,116/2],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',.1,'EdgeAlpha',0);
    %     hold on;line([62,62],get(gca,'YLim'),'color','c')
    %     hold on;line([116/2,116/2],get(gca,'YLim'),'color','k') % temp
    
    title(num2str(NeuronList(NeuronNum)))
    xlim([50 100])
    %ylim([0,1])
end
if includeslc
    legend({'model prediction - control','model prediction - slc', 'true trace - control','true trace - slc'},'Location','northwest');
else
    legend({'model prediction - control', 'true trace - control'},'Location','northwest');
end
%% 0 - 1 - plotting per trial data and mpdel prediction - control only

%%% params
animali = 1;
leaveout = 'neuron'; % neuron or area: all V1
CntrlOnly = 0; % 0 or not
%%% load the lono/loao file here
filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
res = load(fullfile(rootdir,filename));
Ypred = res.Ypred;
TrYcntrl = res.TrYcntrl;

figure;
trN = 1;
NeuronList = [5 6 11 14 15 18 19 20 23 25 26];%[5 11 12 15 18 19 20 23 24 25];% 34 for animal 6, 22 qnd 24 for animal 3
for NeuronNum =1%:10
    % trial averaged
    for trN = 1:10
    subplot(10,1,trN);
    
    plot(Ypred(trN+10,:,NeuronList(NeuronNum)),'color',[0.7 0 0])
    hold on;plot(TrYcntrl(trN+10,:,NeuronList(NeuronNum)),'color',[0 0 0])
    hold on;plot(nanmean(TrYcntrl(:,:,NeuronList(NeuronNum)),1),'color',[.7 .7 .7])
    ylim([0,5])
    yl = get(gca,'YLim');
    hold on;patch([116/2,3*116/4,3*116/4,116/2],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',.1,'EdgeAlpha',0);
    
    title(num2str(NeuronList(NeuronNum)))
    xlim([50 100])
    
    end
end
    legend({'model prediction', 'true trace', 'true trace trial av'},'Location','northwest');

%% 1 - LONO Likelihoods - (mostly control trials relevant here)

rng(1,'twister');
CntrlOnly = 0; % 0 or not
timewindowtype = 'vis'; % vis (55:100) or win (through the silencing timewindow, around 150 ms)
% if vis, cntrl is the same for all lags, if win, even cntrl is lag dependent
nfold = 5;

ndims_pca = 16;
showplot_pca = 0;
animalIds = [1 3 4 5 6];
condi = 1; % nan: no condition 1: only relevant for control trials. uses same condition for both av from dara and from prediction

All_pll_slc_fromAvctrl = cell(1,max(animalIds)); % use a randon control trial to predict it, average pll for all control trials
All_pll_slc_fromModelPred = cell(1,max(animalIds)); % same trial prediction
All_pll_slc_fromDataAv = cell(1,max(animalIds)); % from average across trials from data (silecing trials)

All_pll_cntrl_fromAvctrl = cell(1,max(animalIds));
All_pll_cntrl_fromModelPred = cell(1,max(animalIds)); % same trial prediction
All_pll_cntrl_fromDataAv = cell(1,max(animalIds)); % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
All_pll_cntrl_fromLowdim = cell(1,max(animalIds));
leaveout = 'neuron'; % neuron or area: all V1 ** DO not change, for area, run the next section
for animali = animalIds
    % load the lono/loao file
    filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
    res = load(fullfile(rootdir,filename));
    
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
        filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
        filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    end
    mdl = load(filelist(end).name);
    if strcmp(leaveout,'neuron')
        predictedNeurons = 1:size(mdl.Xval.test_seq{1}(1).y,1);
        npneurons = size(mdl.Xval.test_seq{1}(1).y,1);
    end
    pll_slc_fromAvctrl = nan(8,npneurons);
    pll_slc_fromModelPred = nan(8,npneurons);
    pll_slc_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromModelPred = nan(8,npneurons);
    pll_cntrl_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromAvctrl = nan(8,npneurons);
    pll_cntrl_fromLowdim = nan(8,npneurons);
    
    %%% lag independent vars:
    Ypred = res.Ypred;
    TrYcntrl = res.TrYcntrl;
    %%% make low dimensional Y
    %%% this pools timepoints points and trials together as many points
    %%% (~20000) in the neuron dimensional space (better to use averages instead? [todo])
    TrYcntrl_red = reduce_dim(TrYcntrl,ndims_pca,showplot_pca);
    %%%
    
    for tSilencingLag = 1:8
        
        Ypred_slc = res.Lag{tSilencingLag}.Ypred_slc;
        TrY = res.Lag{tSilencingLag}.TrY;
        
        
        if strcmp(timewindowtype,'win')
            if ~CntrlOnly
                timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            else
                timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            end
            timewindow = timewindow(1:end);
        elseif strcmp(timewindowtype,'vis')
            timewindow = 55:100;%timewindow(1);%:116; % the whole time window
        end
        
        %%% PLLs: look into the function
        [pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
            pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim] = ...
            makePLLs(predictedNeurons,tSilencingLag,timewindow,TrY,TrYcntrl,TrYcntrl_red,Ypred_slc,Ypred,condi,...
            pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
            pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim);
        
        
    end
    All_pll_slc_fromAvctrl{animali} =  pll_slc_fromAvctrl;
    All_pll_slc_fromModelPred{animali} = pll_slc_fromModelPred;
    All_pll_slc_fromDataAv{animali} = pll_slc_fromDataAv;
    
    All_pll_cntrl_fromAvctrl{animali} = pll_cntrl_fromAvctrl;
    All_pll_cntrl_fromModelPred{animali} = pll_cntrl_fromModelPred; % same trial prediction
    All_pll_cntrl_fromDataAv{animali} = pll_cntrl_fromDataAv; % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
    All_pll_cntrl_fromLowdim{animali}= pll_cntrl_fromLowdim;
end
%%% plot PLLs per animal exploratory
if 0
    %%% exploratory, per animal
    figure;
    s = subplot(2,2,1);
    plot(pll_cntrl_fromModelPred-pll_cntrl_fromAvctrl,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromAvctrl,2),'k');
    s.Title.String = 'control trials from model predicion, vs average-control-trial model prediction';
    
    s = subplot(2,2,3);
    plot(pll_slc_fromModelPred-pll_slc_fromAvctrl,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromAvctrl,2),'k');
    s.Title.String = 'slc trials from model predicion, vs average-control-trial model prediction';
    
    %%%% How well again data: (probably only positive for lono, not loao and model performnce (not interaction))
    
    s = subplot(2,2,2);
    plot(pll_cntrl_fromModelPred-pll_cntrl_fromDataAv,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromDataAv,2),'k');
    s.Title.String = 'control trials from model predicion, vs from trial-average data estimate';
    
    s = subplot(2,2,4);
    plot(pll_slc_fromModelPred-pll_slc_fromDataAv,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromDataAv,2),'k');
    s.Title.String = 'slc trials from model predicion, vs from trial-average data estimate';
    
    figure;s=subplot(1,1,1);
    plot(pll_cntrl_fromModelPred-pll_cntrl_fromLowdim,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromLowdim,2),'k');
    s.Title.String = 'control trials from model predicion, vs from low-dim trial-average data estimate';
else
    % assumin 'win', same all lags;
    % only for control animals
    perneuron = cellfun(@(x,y) nanmean(x,1) - nanmean(y,1),All_pll_cntrl_fromModelPred(animalIds),All_pll_cntrl_fromDataAv(animalIds),'UniformOutput',0);
    figure;
    for i =1:numel(perneuron)
        % hold on; scatter(i+zeros(size(perneuron{i})),perneuron{i})
        s = subplot(1,5,i);boxplot(perneuron{i})
        text(1,0.05,['sign-rank p = ',num2str(signrank(perneuron{i}))])
        ylim([-0.05,0.08]);
        s.Title.String = 'pll(model) - pll(Avdata)';
    end
    
    perneuron = cellfun(@(x,y)  nanmean(x,1) - nanmean(y,1),All_pll_cntrl_fromModelPred(animalIds),All_pll_cntrl_fromLowdim(animalIds),'UniformOutput',0);
    figure;
    for i =1:numel(perneuron)
        % hold on; scatter(i+zeros(size(perneuron{i})),perneuron{i})
        s = subplot(1,5,i);boxplot(perneuron{i})
        text(1,0.05,['sign-rank p = ',num2str(signrank(perneuron{i}))])
        ylim([-0.05,0.08])
        s.Title.String = 'pll(model) - pll(AvLowdimData)';
    end
    
    
    perneuron = cellfun(@(x,y)  nanmean(x,1) - nanmean(y,1),All_pll_cntrl_fromModelPred(animalIds),All_pll_cntrl_fromAvctrl(animalIds),'UniformOutput',0);
    figure;
    for i =1:numel(perneuron)
        % hold on; scatter(i+zeros(size(perneuron{i})),perneuron{i})
        s = subplot(1,5,i);boxplot(perneuron{i})
        text(1,0.05,['sign-rank p = ',num2str(signrank(perneuron{i}))])
        ylim([-0.05,0.08])
        s.Title.String = 'pll(model) - pll(modelElseTrialsAverage)';
    end
end
%% 2- LOAO: PLL and ds

rng(1,'twister');
CntrlOnly = 0; % 0 or not
timewindowtype = 'vis'; % vis (55:100) or win (through the silencing timewindow, around 150 ms) -- [fix, for d, it is should be win, for pll vis to compare with above(?)]
% if vis, cntrl is the same for all lags, if win, even cntrl is lag dependent
nfold = 5;

ndims_pca = 16;
showplot_pca = 0;
animalIds = [1 3 4 5 6];
condi = 1; % nan: no condition 1: relevant for both control and slencing trials. uses same condition for both av from dara and from prediction -- used in PLL
calc_pll = 0; % if 0, skips pll calculation, and does only d (to save time)

dMethod = 'mean'; % mean or lda
activity_thresh = nan;% nan or e.g. 0.3750
nboot = 100;%100
calcXd = 1; % 0 or 1, if 1, runs makeDsX instead of makeDs

All_pll_slc_fromAvctrl = cell(1,max(animalIds)); % use a randon control trial to predict it, average pll for all control trials
All_pll_slc_fromModelPred = cell(1,max(animalIds)); % same trial prediction
All_pll_slc_fromDataAv = cell(1,max(animalIds)); % from average across trials from data (silecing trials)
All_pll_cntrl_fromAvctrl = cell(1,max(animalIds));
All_pll_cntrl_fromModelPred = cell(1,max(animalIds)); % same trial prediction
All_pll_cntrl_fromDataAv = cell(1,max(animalIds)); % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
All_pll_cntrl_fromLowdim = cell(1,max(animalIds));
All_dhat = cell(1,max(animalIds));
All_dtrue = cell(1,max(animalIds));
All_dhat_boot = cell(1,max(animalIds));
All_dtrue_boot = cell(1,max(animalIds));
All_dhat_chance = cell(1,max(animalIds));
All_dtrue_chance = cell(1,max(animalIds));
leaveout = 'area'; % neuron or area: all V1 ** Do not change, for neuron, run previous section
for animali = animalIds
    % load the lono/loao file
    filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
    res = load(fullfile(rootdir,filename));
    
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
        filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
        filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    end
    mdl = load(filelist(end).name);
    
    if strcmp(leaveout,'area')
        NV1Cells = res.NV1Cells;
        predictedNeurons = 1:NV1Cells;
        npneurons = NV1Cells;
    end
    
    clear dhat dtrue
    pll_slc_fromAvctrl = nan(8,npneurons);
    pll_slc_fromModelPred = nan(8,npneurons);
    pll_slc_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromModelPred = nan(8,npneurons);
    pll_cntrl_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromAvctrl = nan(8,npneurons);
    pll_cntrl_fromLowdim = nan(8,npneurons);
    if ~calcXd
        dhat = nan(8,npneurons);
        dtrue = nan(8,npneurons);
    else
        dhat.p1 = nan(8,npneurons,nboot);
        dtrue.p1 = nan(8,npneurons,nboot);
        dhat.p2= nan(8,npneurons,nboot);
        dtrue.p2 = nan(8,npneurons,nboot);
    end    
    dhat_boot = nan(8,npneurons,nboot);
    dtrue_boot = nan(8,npneurons,nboot);
    dhat_chance = nan(8,npneurons,nboot);
    dtrue_chance = nan(8,npneurons,nboot);
    %%% lag independent vars:
    Ypred = res.Ypred;
    TrYcntrl = res.TrYcntrl;
    %%% make low dimensional Y
    %%% this pools timepoints points and trials together as many points
    %%% (~20000) in the neuron dimensional space (better to use averages instead? [todo])
    TrYcntrl_red = reduce_dim(TrYcntrl,ndims_pca,showplot_pca);
    %%%
    
    for tSilencingLag = 1:8
        
        Ypred_slc = res.Lag{tSilencingLag}.Ypred_slc;
        TrY = res.Lag{tSilencingLag}.TrY;
        
        
        if strcmp(timewindowtype,'win')
            if ~CntrlOnly
                timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            else
                timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            end
            timewindow = timewindow(1:end);
        elseif strcmp(timewindowtype,'vis')
            timewindow = 55:100;%timewindow(1);%:116; % the whole time window
        end
        
        if strcmp(leaveout, 'area') % ** for calculating d, only 'Win' timewindow is used, timewindowtype is for pll only
            if ~CntrlOnly
                timewindow_d = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            else
                timewindow_d = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
            end
            timewindow_d = timewindow_d(1:end); % end or end-1
            %%% *** put condition like pll, that the trilas without
            %%% prediction are also removed from trials -- there are no
            %%% trials without prediction in loao animals 1 3 4 5 6
            if ~ calcXd
                [dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot] = ...
                    makeDs(dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot,nboot,timewindow_d,NV1Cells,Ypred,Ypred_slc,TrYcntrl,TrY,dMethod,activity_thresh,tSilencingLag);
            else
                [dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot] = ...
                    makeDsX(dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot,nboot,timewindow_d,NV1Cells,Ypred,Ypred_slc,TrYcntrl,TrY,dMethod,activity_thresh,tSilencingLag);               
            end
        end
        
        %%% PLLs
        if calc_pll
            [pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
                pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim] = ...
                makePLLs(predictedNeurons,tSilencingLag,timewindow,TrY,TrYcntrl,TrYcntrl_red,Ypred_slc,Ypred,condi,...
                pll_slc_fromAvctrl,pll_slc_fromModelPred,pll_slc_fromDataAv,...
                pll_cntrl_fromAvctrl,pll_cntrl_fromModelPred,pll_cntrl_fromDataAv,pll_cntrl_fromLowdim);
        end
        
    end
    All_pll_slc_fromAvctrl{animali} =  pll_slc_fromAvctrl;
    All_pll_slc_fromModelPred{animali} = pll_slc_fromModelPred;
    All_pll_slc_fromDataAv{animali} = pll_slc_fromDataAv;
    
    All_pll_cntrl_fromAvctrl{animali} = pll_cntrl_fromAvctrl;
    All_pll_cntrl_fromModelPred{animali} = pll_cntrl_fromModelPred; % same trial prediction
    All_pll_cntrl_fromDataAv{animali} = pll_cntrl_fromDataAv; % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
    All_pll_cntrl_fromLowdim{animali}= pll_cntrl_fromLowdim;
    
    All_dhat{animali} = dhat;
    All_dtrue{animali} = dtrue;
    All_dhat_boot{animali} = dhat_boot;
    All_dtrue_boot{animali} = dtrue_boot;
    All_dhat_chance{animali} = dhat_chance;
    All_dtrue_chance{animali} = dtrue_chance;
end

%%% remove nonexistent animals
All_dhat = All_dhat(animalIds);
All_dtrue = All_dtrue(animalIds);
All_dhat_boot = All_dhat_boot(animalIds);
All_dtrue_boot = All_dtrue_boot(animalIds);
All_dhat_chance = All_dhat_chance(animalIds);
All_dtrue_chance = All_dtrue_chance(animalIds);

if 0
    %%% exploratory, per animal
    figure;
    s = subplot(2,2,1);
    plot(pll_cntrl_fromModelPred-pll_cntrl_fromAvctrl,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromAvctrl,2),'k');
    s.Title.String = 'control trials from model predicion, vs average-control-trial model prediction';
    
    s = subplot(2,2,3);
    plot(pll_slc_fromModelPred-pll_slc_fromAvctrl,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromAvctrl,2),'k');
    s.Title.String = 'slc trials from model predicion, vs average-control-trial model prediction';
    
    %%%% How well again data: (probably only positive for lono, not loao and model performnce (not interaction))
    
    s = subplot(2,2,2);
    plot(pll_cntrl_fromModelPred-pll_cntrl_fromDataAv,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromDataAv,2),'k');
    s.Title.String = 'control trials from model predicion, vs from trial-average data estimate';
    
    s = subplot(2,2,4);
    plot(pll_slc_fromModelPred-pll_slc_fromDataAv,'color',[.8 .8 .8]);
    hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromDataAv,2),'k');
    s.Title.String = 'slc trials from model predicion, vs from trial-average data estimate';
    
    % figure;s=subplot(1,1,1);
    % plot(pll_cntrl_fromModelPred-pll_cntrl_fromLowdim,'color',[.8 .8 .8]);
    % hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromLowdim,2),'k');
    % s.Title.String = 'control trials from model predicion, vs from low-dim trial-average data estimate';
else
    if calc_pll
        perneuron = cellfun(@(x,y) x - y,All_pll_slc_fromModelPred(animalIds),All_pll_slc_fromAvctrl(animalIds),'UniformOutput',0);
        figure;
        for i =1:numel(perneuron)
            % hold on; scatter(i+zeros(size(perneuron{i})),perneuron{i})
            s = subplot(numel(perneuron),1,i);boxplot(perneuron{i}')
            for l = 1:8
                text(l,0.07,num2str(signrank(perneuron{i}(l,:))))
            end
            ylim([-0.05,0.08]);
            hold on;line([0 9], [0 0],'Color','k','LineStyle','--')
            s.Title.String = 'pll(model) - pll(TrialAvCntrlModdel) - pvals signrank';
        end
    end
end

%% d plots
if ~calcXd
    %%% normalize
    for a = 1:numel(All_dhat)
        for i = 1:8
            All_dhat{a}(i,:) = All_dhat{a}(i,:)/norm(All_dhat{a}(i,:));
            All_dtrue{a}(i,:) = All_dtrue{a}(i,:)/norm(All_dtrue{a}(i,:));
            for b = 1:nboot
                All_dhat_boot{a}(i,:,b) = All_dhat_boot{a}(i,:,b)/norm(squeeze(All_dhat_boot{a}(i,:,b)));
                All_dtrue_boot{a}(i,:,b) = All_dtrue_boot{a}(i,:,b)/norm(squeeze(All_dtrue_boot{a}(i,:,b)));
                
                All_dhat_chance{a}(i,:,b) = All_dhat_chance{a}(i,:,b)/norm(squeeze(All_dhat_chance{a}(i,:,b)));
                All_dtrue_chance{a}(i,:,b) = All_dtrue_chance{a}(i,:,b)/norm(squeeze(All_dtrue_chance{a}(i,:,b)));
            end
        end
    end
    %%% calculate angles - compare with predictions from other lag
    
    figure;
    self = [];others = [];
    for a = 1:numel(All_dhat)
        subplot(numel(All_dhat),1,a);
        ch = nan(8,nboot);
        for i = 1:8
            self(end+1) = All_dtrue{a}(i,:)*All_dhat{a}(i,:)';
            hold on; scatter(i,All_dtrue{a}(i,:)*All_dhat{a}(i,:)',100,'k.')
            for j = 1:8
                if j ~= i & abs(j-i)>1
                    others(end+1) = All_dtrue{a}(i,:)*All_dhat{a}(j,:)';
                    hold on; scatter(i+0.1,All_dtrue{a}(i,:)*All_dhat{a}(j,:)','r.')
                end
            end
            for b = 1:nboot
                ch(i,b) = All_dtrue{a}(i,:)*All_dtrue_chance{a}(i,:,b)';
            end
        end
        hold on;line(1:8,prctile(ch',97.5),'Color','k','LineStyle','--');
        hold on;line(1:8,prctile(ch',2.5),'Color','k','LineStyle','--');
    end
    
    figure; histogram(others,[-1:.3:1],'Normalization','pdf')
    hold on; histogram(self,[-1:.3:1],'Normalization','pdf')
    ranksum(self,others)
    %%% calculate angles - compare with data from otherlags
    
    figure;
    for a = 1:numel(All_dhat)
        subplot(numel(All_dhat),1,a);
        for i = 1:8
            hold on; scatter(i,All_dtrue{a}(i,:)*All_dhat{a}(i,:)',100,'k.')
            for j = 1:8
                if j ~= i & abs(j-i)>1
                    hold on; scatter(i+0.1,All_dtrue{a}(i,:)*All_dtrue{a}(j,:)','g.')
                end
            end
        end
        
    end
    %%% show angle decay over time
    figure;
    for a = 1:numel(All_dhat)
        subplot(numel(All_dhat),2,2*a-1);
        mat = All_dtrue{a}*All_dtrue{a}';
        for i = 0:7
            hold on;scatter(i+zeros(size(diag(mat,i))),diag(mat,i),'k.')
        end
        
        subplot(numel(All_dhat),2,2*a);
        mat = All_dhat{a}*All_dhat{a}';
        for i = 0:7
            hold on;scatter(i+zeros(size(diag(mat,i))),diag(mat,i),'c.')
        end
    end
else %%% X case
    %% Normalize
    for a = 1:numel(All_dhat)
        for i = 1:8
            for b = 1:nboot
                All_dhat{a}.p1(i,:,b) = All_dhat{a}.p1(i,:,b)/norm(squeeze(All_dhat{a}.p1(i,:,b)));
                All_dtrue{a}.p1(i,:,b) = All_dtrue{a}.p1(i,:,b)/norm(squeeze(All_dtrue{a}.p1(i,:,b)));
                All_dhat{a}.p2(i,:,b) = All_dhat{a}.p2(i,:,b)/norm(squeeze(All_dhat{a}.p2(i,:,b)));
                All_dtrue{a}.p2(i,:,b) = All_dtrue{a}.p2(i,:,b)/norm(squeeze(All_dtrue{a}.p2(i,:,b)));
                
                
                All_dtrue_chance{a}(i,:,b) = All_dtrue_chance{a}(i,:,b)/norm(squeeze(All_dtrue_chance{a}(i,:,b)));
                All_dhat_chance{a}(i,:,b) = All_dhat_chance{a}(i,:,b)/norm(squeeze(All_dhat_chance{a}(i,:,b)));
            end
        end
    end
    %% self vs others vs chance angles
    %%% 1 - Others consists of lags that are at least minDist lags away from
    %%% the plotte lag.
    %%% ** 2 - Chance level: dhat_chance is a vector going from one half of
    %%% control trials (half of trials) to the other half (half1 assumed to
    %%% be slc, half 2 cntrl)
    %%% Important assumptions for chance: 
    %%% I)the angle is between dtrue and dhat_chance
    %%% II) similar to the non-chance case, p1 (using half 1 cntrl) is
    %%% compared to chane (half 2 cntrl) and p2 is compared to -chance
    %%% (going in opposite direction).
    %%% that's why averahe is negative, not zero
    %%% III) the chance angles are averages over partitioning for each lag (similar to
    %%% other hists) as opposed to distribution across boot (this can be changed by removing mean in line 578)
    minDist = 1;
    figure;
    self = [];others = [];chance = [];
    for a = 1:numel(All_dhat)
        subplot(numel(All_dhat),1,a)
        Xang_ch = nan(8,nboot);
        for i =1:8
            Xang_self = [];
            for b = 1:nboot
                Xang_self(end+1) = (All_dtrue{a}.p1(i,:,b)*All_dhat{a}.p2(i,:,b)' + All_dtrue{a}.p2(i,:,b)*All_dhat{a}.p1(i,:,b)')/2;
            end
            hold on; scatter(i,mean(Xang_self),'k+'); xlim([.5 8.5]);ylim([-1 1]) % average over nboot partitionings is plotted
            self(end+1) = mean(Xang_self);
            
            for j = 1:8
                Xang_cross = [];
                if j ~= i & abs(j-i) > minDist  % lags that are at least minDist away                 
                    for b = 1:nboot
                        Xang_cross(end+1) = (All_dtrue{a}.p1(i,:,b)*All_dhat{a}.p2(j,:,b)' + All_dtrue{a}.p2(i,:,b)*All_dhat{a}.p1(j,:,b)')/2;
                    end
                end
                hold on; scatter(i,mean(Xang_cross),'r.') % average over nboot partitionings is plotted
                others(end+1) = mean(Xang_cross);               
            end
            
            Xang_ch = [];
            for b = 1:nboot
                % Xang_ch(end+1) = (All_dtrue{a}.p1(i,:,b)*All_dtrue_chance{a}(i,:,b)' +  All_dtrue{a}.p2(i,:,b)*(-All_dtrue_chance{a}(i,:,b))')/2;    
                Xang_ch(end+1) = (All_dtrue{a}.p1(i,:,b)*All_dhat_chance{a}(i,:,b)' +  All_dtrue{a}.p2(i,:,b)*(-All_dhat_chance{a}(i,:,b))')/2;   
                % which is the same as : All_dtrue{a}.p1(i,:,b)*All_dhat_chance{a}(i,:,b)';          
                % hold on; scatter(i-0.1,Xang_ch(end),'g.')                
            end           
            chance = [chance mean(Xang_ch)]; % average over nboot partitionings is added             
        end
    end
    
    figure; histogram(others,[-1:.25:1],'Normalization','pdf','FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0);hold on;%line([nanmedian(others),nanmedian(others)],[0 1.5],'Color','r')
    hold on; histogram(self,[-1:.25:1],'Normalization','pdf','FaceColor','k','FaceAlpha',0.5,'EdgeAlpha',0);hold on;%line([nanmedian(self),nanmedian(self)],[0 1.5],'Color','k')
%    hold on; histogram(chance,[-1:.05:1],'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.7 .7 .7])
    hold on; errorbar(prctile(chance,50),0.1,(prctile(chance,50) - prctile(chance,2.5)),(prctile(chance,97.5)-prctile(chance,50)),'horizontal','Color','k','LineWidth',2)
    legend({'self',['others at distance ',num2str(minDist)],'chance'})
    hold on; text(-1,1,['self vs. others ranksum p = ',num2str(ranksum(self,others))])
    hold on; text(-1,0.9,['self median = ',num2str(nanmedian(self))]) 
    hold on; text(-1,0.8,['others median = ',num2str(nanmedian(others))]) 
    
    %% decay over time true vs model prediction
    alltrue =cell(1,8);
    allhat = cell(1,8);
    
    figure;
    for a = 1:numel(All_dhat)
        subplot(numel(All_dhat),2,2*a-1);
        mat = [];
        for b = 1:nboot
            mat = cat(3,mat,All_dtrue{a}.p1(:,:,b)*All_dtrue{a}.p2(:,:,b)');
        end
        mat = nanmean(mat,3);
        x=[];
        g = [];
        for i = 0:7
            hold on;scatter(i+zeros(size(diag(mat,i))),diag(mat,i),'k.')
            x = [x;diag(mat,i)];
            g = [g;i+zeros(size(diag(mat,i)))];
            alltrue{i+1} =  [alltrue{i+1};diag(mat,i)];
        end
        xlim([-0.5,8])
        % hold on;boxplot(x,g)
        
        subplot(numel(All_dhat),2,2*a);
        mat = [];
        for b = 1:nboot
            mat = cat(3,mat,All_dhat{a}.p1(:,:,b)*All_dhat{a}.p2(:,:,b)');
        end
        mat = nanmean(mat,3);
        x=[];
        g = [];
        for i = 0:7
            hold on;scatter(i+zeros(size(diag(mat,i))),diag(mat,i),'b.')
            x = [x;diag(mat,i)];
            g = [g;i+zeros(size(diag(mat,i)))];
            allhat{i+1} =  [allhat{i+1};diag(mat,i)];
        end
        xlim([-0.5,8])
        %  hold on;boxplot(x,g)
    end
    figure
    for i = 1:8
        subplot(1,2,1);hold on; scatter(i-1+zeros(numel(alltrue{i}),1),alltrue{i},'k')
        subplot(1,2,2);hold on; scatter(i-1+zeros(numel(allhat{i}),1),allhat{i},'b')
    end
    
    figure
    for i = 1:8
        subplot(1,2,1);hold on; scatter(i-1,nanmean(alltrue{i}),'k.');hold on; errorbar(i-1,nanmean(alltrue{i}),2*nanstd(alltrue{i})/sqrt(numel(alltrue{i})),'Color','k')
        ylim([-0.7,1]);xlim([-0.5,7.5])
        subplot(1,2,2);hold on; scatter(i-1,nanmean(allhat{i}),'b.');hold on; errorbar(i-1,nanmean(allhat{i}),2*nanstd(allhat{i})/sqrt(numel(allhat{i})),'Color','b')
        ylim([-0.7,1]);xlim([-0.5,7.5])
    end
    subplot(1,2,1);hold on;plot(0:7,cellfun(@(x) nanmean(x),alltrue),'k')
    subplot(1,2,2);hold on;plot(0:7,cellfun(@(x) nanmean(x),allhat),'b')
end