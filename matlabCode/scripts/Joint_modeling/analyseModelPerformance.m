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

%% 1 - LONO Likelihoods - (mostly control trials relevant here)

rng(1,'twister');
CntrlOnly = 0; % 0 or not
leaveout = 'neuron'; % neuron or area: all V1
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
    
    predictedNeurons = 1:size(mdl.Xval.test_seq{1}(1).y,1);
    npneurons = size(mdl.Xval.test_seq{1}(1).y,1);
    
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
        
        %%% PLLs
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
%% LOAO -- haven't started yet -- combine draft_analye_...2 (relevant part already copied here) and the above section
rng(1,'twister');
CntrlOnly = 0; % 0 or not
nboot = 1;
nfold = 5;

leaveout = 'area'; % neuron or area: all V1
dMethod = 'mean'; % mean or lda
activity_thresh = nan;% nan or e.g. 0.3750

if CntrlOnly
    UInd = 1;
else
    UInd = 1:9;
end
for animali = 1:6
    % load the lono/loao file
    
    if strcmp(leaveout, 'area')
        predictedNeurons = 1:NV1Cells;
        npneurons = NV1Cells;
    else
        predictedNeurons = 1:size(res.Xval.test_seq{fold}(1).y,1);
        npneurons = size(res.Xval.test_seq{fold}(1).y,1);
    end
    
    
    dhat = nan(8,npneurons);
    dtrue = nan(8,npneurons);
    dhat_boot = nan(8,npneurons,nboot);
    dtrue_boot = nan(8,npneurons,nboot);
    dhat_chance = nan(8,npneurons,nboot);
    dtrue_chance = nan(8,npneurons,nboot);
    pll_slc_fromAvctrl = nan(8,npneurons);
    pll_slc_fromModelPred = nan(8,npneurons);
    pll_slc_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromModelPred = nan(8,npneurons);
    pll_cntrl_fromDataAv = nan(8,npneurons);
    pll_cntrl_fromAvctrl = nan(8,npneurons);
    
    for tSilencingLag = 1:8
        
        if ~CntrlOnly
            timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
        else
            timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, res.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
        end
        % optional: shorted window
        timewindow = timewindow(1:end);
        
        if strcmp(leaveout, 'area')
            [dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot] = ...
                makeDs(dhat,dtrue,dhat_chance,dtrue_chance,dhat_boot,dtrue_boot,nboot,timewindow,NV1Cells,Ypred,Ypred_slc,TrYcntrl,TrY,dMethod,activity_thresh,tSilencingLag);
        end
        % ooption to use a different time window for pll
        timewindow = 55:100;%timewindow(1);%:116; % the whole time window
        
        
        
        %%% pll
        for NeuronNum = predictedNeurons
            
            pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = 0; % use a randon control trial to predict it, average pll for all control trials
            pll_slc_fromModelPred(tSilencingLag,NeuronNum) = 0; % same trial prediction
            pll_slc_fromDataAv(tSilencingLag,NeuronNum) = 0; % from average across trials from data (silecing trials)
            
            pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = 0;
            pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = 0; % same trial prediction
            pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) = 0; % from average across trials from data (control trials) -- if 0 returns inf and nan, fix [TODO]
            % also, the trial itseld should be excluded, average others  (for both model control and av data) [TODO]
            
            
            for slctr = 1:size(TrY(:,timewindow,NeuronNum),1)
                elsetr = setdiff(1:size(TrY(:,timewindow,NeuronNum),1),slctr);
                for timepoint = 1:numel(timewindow)
                    if Ypred_slc(slctr,timewindow(timepoint),NeuronNum)>0
                        pll_slc_fromModelPred(tSilencingLag,NeuronNum) = nansum([pll_slc_fromModelPred(tSilencingLag,NeuronNum) , ...
                            TrY(slctr,timewindow(timepoint),NeuronNum)*log(Ypred_slc(slctr,timewindow(timepoint),NeuronNum))-Ypred_slc(slctr,timewindow(timepoint),NeuronNum)]);
                    end
                    if Ypred(:,timewindow(timepoint),NeuronNum)>0
                        pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = nansum([pll_slc_fromAvctrl(tSilencingLag,NeuronNum) , ...
                            mean(TrY(slctr,timewindow(timepoint),NeuronNum).*log(Ypred(:,timewindow(timepoint),NeuronNum))-Ypred(:,timewindow(timepoint),NeuronNum))]);
                    end
                    if mean(TrY(elsetr,timewindow(timepoint),NeuronNum),1)>0
                        pll_slc_fromDataAv(tSilencingLag,NeuronNum) = nansum([pll_slc_fromDataAv(tSilencingLag,NeuronNum) , ...
                            TrY(slctr,timewindow(timepoint),NeuronNum).*log(mean(TrY(elsetr,timewindow(timepoint),NeuronNum),1))-mean(TrY(elsetr,timewindow(timepoint),NeuronNum),1)]);
                    end
                end
            end
            
            % average per number of trials and time points, to make
            % comparable
            
            pll_slc_fromModelPred(tSilencingLag,NeuronNum) = pll_slc_fromModelPred(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
            pll_slc_fromAvctrl(tSilencingLag,NeuronNum) = pll_slc_fromAvctrl(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
            pll_slc_fromDataAv(tSilencingLag,NeuronNum) = pll_slc_fromDataAv(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrY(:,timewindow,NeuronNum),1));
            
            for ctrtr = 1:size(TrYcntrl(:,timewindow,NeuronNum),1)
                elsetr = setdiff(1:size(TrYcntrl(:,timewindow,NeuronNum),1),ctrtr);
                for timepoint = 1:numel(timewindow)
                    if Ypred(ctrtr,timewindow(timepoint),NeuronNum)>0
                        pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) , ...
                            TrYcntrl(ctrtr,timewindow(timepoint),NeuronNum)*log(Ypred(ctrtr,timewindow(timepoint),NeuronNum))-Ypred(ctrtr,timewindow(timepoint),NeuronNum)]);
                    end
                    if Ypred(elsetr,timewindow(timepoint),NeuronNum)>0
                        pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) , ...
                            mean(TrYcntrl(ctrtr,timewindow(timepoint),NeuronNum).*log(Ypred(elsetr,timewindow(timepoint),NeuronNum))-Ypred(elsetr,timewindow(timepoint),NeuronNum))]);
                    end
                    if mean(TrYcntrl(elsetr,timewindow(timepoint),NeuronNum),1)>0
                        pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) , ...
                            TrYcntrl(ctrtr,timewindow(timepoint),NeuronNum).*log(mean(TrYcntrl(elsetr,timewindow(timepoint),NeuronNum),1))-mean(TrYcntrl(elsetr,timewindow(timepoint),NeuronNum),1)]);
                    end
                    %                     % pca
                    %                     if mean(TrYcntrl_red(elsetr,timewindow(timepoint),NeuronNum),1)>0
                    %                         pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) = nansum([pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) , ...
                    %                             TrYcntrl(ctrtr,timewindow(timepoint),NeuronNum).*log(mean(TrYcntrl_red(elsetr,timewindow(timepoint),NeuronNum),1))-mean(TrYcntrl_red(elsetr,timewindow(timepoint),NeuronNum),1)]);
                    %                     end
                end
            end
            pll_cntrl_fromModelPred(tSilencingLag,NeuronNum) = pll_cntrl_fromModelPred(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
            pll_cntrl_fromDataAv(tSilencingLag,NeuronNum) =  pll_cntrl_fromDataAv(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
            pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum) = pll_cntrl_fromAvctrl(tSilencingLag,NeuronNum)/(numel(timewindow)*size(TrYcntrl(:,timewindow,NeuronNum),1));
        end
        
    end
    
end
%% PLLs
figure;
s = subplot(2,2,1);
plot(pll_cntrl_fromModelPred(2,:)-pll_cntrl_fromAvctrl(2,:),'color',[.8 .8 .8]);
hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromAvctrl,2),'k');
s.Title.String = 'control trials from model predicion, vs average-control-trial model prediction';

s = subplot(2,2,3);
plot(pll_slc_fromModelPred(2,:)-pll_slc_fromAvctrl(2,:),'color',[.8 .8 .8]);
hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromAvctrl,2),'k');
s.Title.String = 'slc trials from model predicion, vs average-control-trial model prediction';

%%%% How well again data: (probably only positive for lono, not loao and model performnce (not interaction))

s = subplot(2,2,2);
plot(pll_cntrl_fromModelPred(2,:)-pll_cntrl_fromDataAv(2,:),'color',[.8 .8 .8]);
hold on;plot(mean(pll_cntrl_fromModelPred-pll_cntrl_fromDataAv,2),'k');
s.Title.String = 'control trials from model predicion, vs from trial-average data estimate';

s = subplot(2,2,4);
plot(pll_slc_fromModelPred(2,:)-pll_slc_fromDataAv(2,:),'color',[.8 .8 .8]);
hold on;plot(mean(pll_slc_fromModelPred-pll_slc_fromDataAv,2),'k');
s.Title.String = 'slc trials from model predicion, vs from trial-average data estimate';
%% PCA scratch
%X = (squeeze(nanmean(TrY(:,timewindow,:),1)));%timepoints*neurons
X = reshape(TrY,[],size(TrY,3));%timepointsandtrials*neurons

[coefs,scores,~] = pca(X,'Centered',0);
% X = scores*coefs'
coefs_red = [coefs(:,1:16),zeros(size(coefs,1),size(coefs,2)-16)]; % assuming descending order, mean subtraction?
Xred = scores*coefs_red';
%figure;imagesc(X);figure;imagesc(Xred)

TrYOut=reshape(Xred,size(TrY,1),size(TrY,2),size(TrY,3));
figure;imagesc(squeeze(nanmean(TrY(:,:,:),1)));figure;imagesc(squeeze(nanmean(TrYOut(:,:,:),1)))

%%%%%
X = reshape(TrYcntrl,[],size(TrY,3));%timepointsandtrials*neurons

[coefs,scores,~] = pca(X,'Centered',0);
% X = scores*coefs'
coefs_red = [coefs(:,1:16),zeros(size(coefs,1),size(coefs,2)-16)]; % assuming descending order, mean subtraction?
Xred = scores*coefs_red';
%figure;imagesc(X);figure;imagesc(Xred)

TrYcntrl_red=reshape(Xred,size(TrYcntrl,1),size(TrYcntrl,2),size(TrYcntrl,3));
figure;imagesc(squeeze(nanmean(TrYcntrl(:,:,:),1)));figure;imagesc(squeeze(nanmean(TrYcntrl_red(:,:,:),1)))
%% before normalizing: distribution of norms
magdist = nan(8,nboot);
mag = nan(1,8);
zs = nan(1,8);
for i = 1:8
    mag(i)=norm(dhat(i,:));
    for b = 1:nboot
        magdist(i,b) = norm(squeeze(dhat_chance(i,:,b)));
    end
    [i,(mag(i) - mean(magdist(i,:)))/std(magdist(i,:))]
    zs(i) = (mag(i) - mean(magdist(i,:)))/std(magdist(i,:));
end
figure;plot(zs)
mag = nan(1,8);
maghat = nan(1,8);
for i = 1:8
    mag(i)=norm(dtrue(i,:));
    maghat(i)=norm(dhat(i,:));
    
end
figure;scatter(mag,maghat)
%% normalize

for i = 1:8
    dhat(i,:) = dhat(i,:)/norm(dhat(i,:));
    dtrue(i,:) = dtrue(i,:)/norm(dtrue(i,:));
    
    for b = 1:nboot
        dhat_boot(i,:,b) = dhat_boot(i,:,b)/norm(squeeze(dhat_boot(i,:,b)));
        dtrue_boot(i,:,b) = dtrue_boot(i,:,b)/norm(squeeze(dtrue_boot(i,:,b)));
        
        dhat_chance(i,:,b) = dhat_chance(i,:,b)/norm(squeeze(dhat_chance(i,:,b)));
    end
end
%% calculate chance level, if time variance was not captured

figure;subplot(1,2,1);
self = [];
others = [];
for i = 1:8
    hold on; scatter(i,dtrue(i,:)*dhat(i,:)',100,'k.')
    self(end+1) = dtrue(i,:)*dhat(i,:)';
    dis = [];
    for j = 1:8
        if j ~= i & abs(j-i)>2
            hold on; scatter(i+0.1,dtrue(i,:)*dhat(j,:)','c.')
            others(end+1) = dtrue(i,:)*dhat(j,:)';
            hold on; scatter(i-0.1,dtrue(i,:)*dtrue(j,:)','g.')
            for b = 1:nboot
                %                  hold on; scatter(i-0.1,dtrue_boot(i,:,b)*dtrue_boot(j,:,b)','g.')
                hold on; scatter(i-0.1,dtrue(i,:)*dtrue_boot(j,:,b)','g.')
                dis(end+1) = dtrue_boot(i,:,b)*dtrue_boot(j,:,b)';
            end
        end
    end
    zs = (dtrue(i,:)*dhat(i,:)' - mean(dis))/std(dis)
end
xlim([0 9]);
subplot(1,2,2);scatter(randn(1,numel(self)),self,100,'k.');
hold on;scatter(10+randn(1,numel(others)),others,100,'c.');
xlim([-5,15])
% 1 and 4 don't work, I think needs noise measure
figure;imagesc(dhat*dhat')
figure;imagesc(dtrue*dtrue')
%% bootstrap
figure;subplot(1,2,1);
boot = nan(8,nboot);
for i = 1:8
    for b = 1:nboot
        hold on; scatter(i,dtrue_boot(i,:,b)*dhat_boot(i,:,b)',100,'k.')
        %hold on; scatter(i+0.1,dtrue_boot(i,:,b)*dtrue_boot(i,:,b)',100,'c.')
        boot(i,b) = dtrue_boot(i,:,b)*dhat_boot(i,:,b)';
    end
    
end
xlim([0 9]);

% 1 and 4 don't work, I think needs noise measure
%figure;imagesc(dhat*dhat')
%figure;imagesc(dtrue*dtrue')
%% chance level
figure;subplot(1,2,1);
ch = nan(8,nboot);
for i = 1:8
    for b = 1:nboot
        ch(i,b) = dtrue(i,:)*dhat_chance(i,:,b)';
        hold on; scatter(i,dtrue(i,:)*dhat_chance(i,:,b)',100,'k.')
    end
end
xlim([0 9]);
%%%
zs = nan(1,8);
for i = 1:8
    zs(i) = (dtrue(i,:)*dhat(i,:)' -  mean(ch(i,:)))/std(ch(i,:))
end
figure;plot(zs)
%% get input timeline
vis = find(sum(cell2mat(cellfun(@(y) y(1,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l1 = find(sum(cell2mat(cellfun(@(y) y(2,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l2 = find(sum(cell2mat(cellfun(@(y) y(3,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l3 = find(sum(cell2mat(cellfun(@(y) y(4,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l4 = find(sum(cell2mat(cellfun(@(y) y(5,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l5 = find(sum(cell2mat(cellfun(@(y) y(6,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l6 = find(sum(cell2mat(cellfun(@(y) y(7,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l7 = find(sum(cell2mat(cellfun(@(y) y(8,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
l8 = find(sum(cell2mat(cellfun(@(y) y(9,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
%%
% only do 1 animal from above now
% direction of influence will be calculatd, instead of LDA, just the
% difference between mus. But it could still have a distribution: half-half
% trials
% for now, just the average
timewindow = l3; % should be same as tSilencing in the above section
timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
dhat = nan(8,27);
dtrue = nan(8,27);
for NeuronNum = 1:27
    dhat(tSilencingLag,NeuronNum) = (nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),1)) - nanmean(nanmean(Ypred(:,timewindow,(NeuronNum)),1)));
    dtrue(tSilencingLag,NeuronNum) = (nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),1)) - nanmean(nanmean(TrYcntrl(:,timewindow,(NeuronNum)),1)));
    %dhat_chance =
end
similarity = (dhat*dtrue')/(norm(dhat)*norm(dtrue))
