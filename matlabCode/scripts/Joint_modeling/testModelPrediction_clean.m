
animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

% orthogonalization not relevant here
% This code can chooose between leave-one-area-out, simialr to analyse_predict_Xtrials
% andleave-one-neuron, as in analyse_predict_Xtrials_LONO, and can
% when slc trials used, both input is received and y of LM is observed
% (which is silenced)

doplot = 0;
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

for animali = 3%:6
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
        filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
        filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    end
    
    
    res = load(filelist(end).name);
    
    NV1Cells = find(sum(res.Xval.params{1}.model.C(:,1:8),2)==0,1)-1%27;% all 37 % make dynamic maybe try animals with more lm cells
    pvcell =[]%36
    
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
        Ypred = cell(1,nfold);
        Ypred_slc = cell(1,nfold);
        TrY = cell(1,nfold);
        TrYcntrl = cell(1,nfold);  
        
        
        % maybe push tSilencinglag inside to reduce 7 extra control
        % calculation, an save the 4 traces afterwards. Then load for the
        % next part.
        for fold = 1:nfold
            if strcmp(leaveout, 'area')
                [Ypred{fold},Ypred_slc{fold},TrYcntrl{fold},TrY{fold}] = makeYs(res,NV1Cells,pvcell,UInd,fold,tSilencingLag);
            elseif strcmp(leaveout, 'neuron')
                warning off % warning of sigular matrix in laplaceinferencecore
                [Ypred{fold},Ypred_slc{fold},TrYcntrl{fold},TrY{fold}] = makeYs_LONO(res,UInd,fold,tSilencingLag);
            end
        end
        %%% here to the number of folds
        
        Ypred = cat(1,Ypred{1},Ypred{2},Ypred{3},Ypred{4},Ypred{5});
        Ypred_slc = cat(1,Ypred_slc{1},Ypred_slc{2},Ypred_slc{3},Ypred_slc{4},Ypred_slc{5});
        TrY = cat(1,TrY{1},TrY{2},TrY{3},TrY{4},TrY{5});
        TrYcntrl = cat(1,TrYcntrl{1},TrYcntrl{2},TrYcntrl{3},TrYcntrl{4},TrYcntrl{5});
        
        %%%% at this point, the real traces and model predictions in X
        %%%% trials and using X neurons are ready, from here on,
        %%%% performance meansures are calculated
        
        if doplot
            figure;
            NeuronList =11:40;%[5 11 12 15 18 19 20 23 24 25];% 34 for animal 6, 22 qnd 24 for animal 3
            for NeuronNum = 1:10%1:27
                % trial averaged
                subplot(10,1,NeuronNum);
                plot(nanmean(Ypred(:,:,NeuronList(NeuronNum)),1),'color',[.7 .7 .7])
                hold on;plot(nanmean(Ypred_slc(:,:,NeuronList(NeuronNum)),1),'color','c')
                hold on;plot(nanmean(TrYcntrl(:,:,NeuronList(NeuronNum)),1),'k')
                hold on;plot(nanmean(TrY(:,:,NeuronList(NeuronNum)),1),'b')
                hold on;line([62,62],get(gca,'YLim'),'color','c')
                hold on;line([116/2,116/2],get(gca,'YLim'),'color','k') % temp
                
                title(num2str(NeuronList(NeuronNum)))
                xlim([50 100])
                %ylim([0,1])
            end
            
            legend({'model prediction - control','model prediction - slc', 'tru trace - control','true trace - slc'},'Location','northwest');
        end
        
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
