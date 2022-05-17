animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

doplot = 0;
rng(1,'twister');
% orthogonalization not relevant here
% This is leave-one-area-out, simialr to analyse_predict_Xtrials
% Later add leave-one-neuron out for comparsion, as in
% analyyse_predict_Xtrials_LONO

CntrlOnly = 0; % 0 or not


if CntrlOnly
    UInd = 1;
else
    UInd = 1:9;
end

nboot = 1;

dhat = nan(8,27);
dtrue = nan(8,27);
dhat_boot = nan(8,27,nboot);
dtrue_boot = nan(8,27,nboot);
dhat_chance = nan(8,27,nboot);
dtrue_chance = nan(8,27,nboot);

pll_ctrl = nan(8,27);
pll_slc = nan(8,27);
mse_ctrl = nan(8,27);
mse_slc = nan(8,27);
% dhat = nan(8,27);
% dtrue = nan(8,27);

for animali  =1%:6
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    end
    
    
    filelist = dir(['*splt1_onlyCorrect_exGrooming_go_5FoldXval*.mat']);
    res = load(filelist(end).name);
    nfold = 5;
    
    for   tSilencingLag = 1:8
        Ypred = cell(1,nfold);
        Ypred_slc = cell(1,nfold);
        TrY = cell(1,nfold);
        TrYcntrl = cell(1,nfold);
        
        for    fold = 1:nfold
            
            NV1Cells = find(sum(res.Xval.params{1}.model.C(:,1:8),2)==0,1)-1%27;% all 37 % make dynamic maybe try animals with more lm cells
            pvcell =[]%36
            % lag = 2;
            
            %%% these all from test trials
            TrU = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
            TrUv = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
            TrY{fold} = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,NV1Cells);
            TrYcntrl{fold} = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,NV1Cells);
            seqLM_slc = [];  % from test trials
            seqLM = [];
            TrNum =[];
            TrNumcntrl = [];
            for tr = 1:numel(res.Xval.test_seq{fold})
                if sum(sum(res.Xval.test_seq{fold}(tr).u(2:end,:)))==0
                    
                    TrYcntrl{fold}(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(1:NV1Cells,:)';
                    TrNumcntrl(end+1) = tr;
                    seqLM(end+1).y = res.Xval.test_seq{fold}(tr).y(NV1Cells+1:end,:);
                    seqLM(end).y(pvcell-NV1Cells,:) = [];
                    seqLM(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
                    seqLM(end).T = res.Xval.test_seq{fold}(tr).T;
                    %  elseif(slcres.seq(tr).u(2,11) - slcres.seq(tr).u(2,10)) > 0 % change to else, if you want all silencing onsets not only lag 2
                    % elseif(res.Xval.test_seq{fold}(tr).u(2,62) - res.Xval.test_seq{fold}(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
                elseif sum(sum(res.Xval.test_seq{fold}(tr).u(1+tSilencingLag,:))) > 0
                    %elseif(res.Xval.test_seq{fold}(tr).u(2,85) - res.Xval.test_seq{fold}(tr).u(2,84)) > 0
                    seqLM_slc(end+1).y = res.Xval.test_seq{fold}(tr).y(NV1Cells+1:end,:);
                    seqLM_slc(end).y(pvcell-NV1Cells,:) = [];
                    seqLM_slc(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
                    %                 seqLM_slc(end).u = [res.Xval.test_seq{fold}(tr).u(1,:);...
                    %                     zeros(size(res.Xval.test_seq{fold}(tr).u(2,:)))]; % 1 0r : if res= slcres
                    seqLM_slc(end).T = res.Xval.test_seq{fold}(tr).T;
                    
                    TrNum(end+1) = tr;
                    TrU(tr,:) =  res.Xval.test_seq{fold}(tr).u(1+tSilencingLag,:); %% check this, modified
                    TrUv(tr,:) =  res.Xval.test_seq{fold}(tr).u(1,:);
                    TrY{fold}(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(1:NV1Cells,:)';
                end
            end
            
            if 0 % why ? check analyse_predict_xtrials
                u = seqLM_slc(1).u(:,:);
                clear seqLM_slc;
                seqLM_slc = seqLM;
                for tr =1:length(seqLM_slc)
                    seqLM_slc(tr).u(:,:) = u;
                end
            end
            
            
            TrU = TrU(TrNum,:);
            TrUv = TrUv(TrNum,:);
            TrY{fold} = TrY{fold}(TrNum,:,:);
            TrYcntrl{fold} = TrYcntrl{fold}(TrNumcntrl,:,:);
            
            paramsLM = res.Xval.params{fold};
            paramsLM.model.C(pvcell,:)=[];
            paramsLM.model.d(pvcell,:)=[];
            paramsLM.model.C =  paramsLM.model.C(NV1Cells+1:end,:);
            paramsLM.model.d=  paramsLM.model.d(NV1Cells+1:end);
            
            [seqLM_slc,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM_slc);
            [seqLM,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM);
            
            Ypred{fold} = nan(length(seqLM),res.Xval.train_seq{fold}(1).T,NV1Cells);
            Ypred_slc{fold} = nan(length(seqLM_slc),res.Xval.train_seq{fold}(1).T,NV1Cells);
            
            C = res.Xval.params{fold}.model.C(1:NV1Cells,1:8);
            d = res.Xval.params{fold}.model.d(1:NV1Cells);
            for tr = 1:length(seqLM_slc)
                if 0
                    Ypred_slc{fold}(tr,:,:) = exp(C*seqLM_slc(tr).posterior.xsm(1:8,:) + d)';
                    %Ypred_slc{fold}(tr,:,:) = (C*seqLM_slc(tr).posterior.xsm(1:8,:) + d)';
                else
                    Vsm   = reshape(seqLM_slc(tr).posterior.Vsm' ,16,16,res.Xval.train_seq{fold}(1).T);
                    for t = 1:res.Xval.train_seq{fold}(1).T
                        Sigma = Vsm(1:8,1:8,t); % are the first 2 correct?
                        Ypred_slc{fold}(tr,t,:) = (exp(C*seqLM_slc(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)');
                    end
                end
            end
            
            for tr = 1:length(seqLM)
                if 0
                    Ypred{fold}(tr,:,:) = exp(C*seqLM(tr).posterior.xsm(1:8,:) + d)';
                    %Ypred{fold}(tr,:,:) = (C*seqLM(tr).posterior.xsm(1:8,:) + d)';
                else
                    Vsm   = reshape(seqLM(tr).posterior.Vsm' ,16,16,res.Xval.train_seq{fold}(1).T);
                    for t = 1:res.Xval.train_seq{fold}(1).T
                        Sigma = Vsm(1:8,1:8,t); % are the first 2 correct?
                        Ypred{fold}(tr,t,:) = (exp(C*seqLM(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)'); % check diag
                    end
                end
                
            end
            
            
        end
        
        %     % here to the number of folds
        %     Ypred = cat(1,Ypred{1},Ypred{2});
        %     Ypred_slc = cat(1,Ypred_slc{1},Ypred_slc{2});
        %     TrY = cat(1,TrY{1},TrY{2});
        %     TrYcntrl = cat(1,TrYcntrl{1},TrYcntrl{2});
        
        Ypred = cat(1,Ypred{1},Ypred{2},Ypred{3},Ypred{4},Ypred{5});
        Ypred_slc = cat(1,Ypred_slc{1},Ypred_slc{2},Ypred_slc{3},Ypred_slc{4},Ypred_slc{5});
        TrY = cat(1,TrY{1},TrY{2},TrY{3},TrY{4},TrY{5});
        TrYcntrl = cat(1,TrYcntrl{1},TrYcntrl{2},TrYcntrl{3},TrYcntrl{4},TrYcntrl{5});
        
        if doplot
            figure;
            NeuronList =[5 11 12 15 18 19 20 23 24 25];% 34 for animal 6, 22 qnd 24 for animal 3
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
        
        timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, res.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
        
        % optional: shorted window
        timewindow = timewindow(1:end);
       % timewindow = 1:timewindow(1);%:116; % the whole time window
        for NeuronNum = 1:NV1Cells
            dhat(tSilencingLag,NeuronNum) = nanmean(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(:,timewindow,(NeuronNum)),2));
            dtrue(tSilencingLag,NeuronNum) = nanmean(nanmean(TrY(:,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(:,timewindow,(NeuronNum)),2));
            
            
            for b = 1:nboot
                
               
                slc_trials = datasample(1:numel(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)),numel(nanmean(Ypred_slc(:,timewindow,(NeuronNum)),2)));
                cntrl_trials = datasample(1:numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)),numel(nanmean(Ypred(:,timewindow,(NeuronNum)),2)));
                
                dhat_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred(cntrl_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(:,timewindow,(NeuronNum)),2));
                dtrue_chance(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrYcntrl(cntrl_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(:,timewindow,(NeuronNum)),2));
                
                dhat_boot(tSilencingLag,NeuronNum,b) = nanmean(nanmean(Ypred_slc(slc_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(Ypred(cntrl_trials,timewindow,(NeuronNum)),2));                
                dtrue_boot(tSilencingLag,NeuronNum,b) = nanmean(nanmean(TrY(slc_trials,timewindow,(NeuronNum)),2)) - nanmean(nanmean(TrYcntrl(cntrl_trials,timewindow,(NeuronNum)),2));
            end
            
            
            %%% mse
            % mse_ctrl(tSilencingLag,NeuronNum) is error between
            % trial-averaged predictions of cntrl trial and the real
            % silencing trial
            
            % nanmean(Ypred(:,timewindow,(NeuronNum)),1) % trial averaged
            
            pll_ctrl(tSilencingLag,NeuronNum) = 0;
            pll_slc(tSilencingLag,NeuronNum) = 0;
%             mse_ctrl(tSilencingLag,NeuronNum) = 0;
%             mse_slc(tSilencingLag,NeuronNum) = 0;
            
               % mse
%             mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum) + mean((nanmean(TrY(:,timewindow,(NeuronNum)),1) - nanmean(Ypred(:,timewindow,(NeuronNum)),1)).^2) ;% this mean is average over time
%             mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum) + mean((nanmean(TrY(:,timewindow,(NeuronNum)),1) - nanmean(Ypred_slc(:,timewindow,(NeuronNum)),1)).^2) ;% this mean is average over time

            for slctr = 1:min(size(TrY(:,timewindow,(NeuronNum)),1),size(TrYcntrl(:,timewindow,(NeuronNum)),1))-1
                % mse
%                  mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum) + mean((TrY(slctr,timewindow,(NeuronNum)) - nanmean(Ypred(:,timewindow,(NeuronNum)),1)).^2) ;% this mean is average over time
%                  mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum) + mean((TrY(slctr,timewindow,(NeuronNum)) - nanmean(Ypred_slc(:,timewindow,(NeuronNum)),1)).^2) ;% this mean is average over time
%                mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum) + mean((TrY(slctr,timewindow,NeuronNum) - nanmean(Ypred(:,timewindow,NeuronNum),1)).^2) ;% this mean is average over time
%                mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum) + mean((TrY(slctr,timewindow,NeuronNum) - Ypred_slc(slctr,timewindow,NeuronNum)).^2) ;% this mean is average over time
%                 mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum) + mean((nanmean(TrY(:,timewindow,NeuronNum),1) - Ypred(slctr,timewindow,NeuronNum)).^2) ;% this mean is average over time
%                 mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum) + mean((nanmean(TrY(:,timewindow,NeuronNum),1) - Ypred_slc(slctr,timewindow,NeuronNum)).^2) ;% this mean is average over time

%             % poisson logprobability
             for timepoint = 1:numel(timewindow)
%                 mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum) + log(poisspdf(TrY(slctr,timewindow(timepoint),(NeuronNum)), nanmean(Ypred(:,timewindow(timepoint),(NeuronNum)),1)));
%                 mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum) + log(poisspdf(TrY(slctr,timewindow(timepoint),(NeuronNum)), nanmean(Ypred_slc(:,timewindow(timepoint),(NeuronNum)),1)));

%                  pll_ctrl(tSilencingLag,NeuronNum) = pll_ctrl(tSilencingLag,NeuronNum) + TrY(slctr,timewindow(timepoint),(NeuronNum))*log(Ypred(slctr,timewindow(timepoint),NeuronNum))-Ypred(slctr,timewindow(timepoint),NeuronNum);
                if 1
                pll_slc(tSilencingLag,NeuronNum) = pll_slc(tSilencingLag,NeuronNum) + TrY(slctr,timewindow(timepoint),(NeuronNum))*log(Ypred_slc(slctr,timewindow(timepoint),NeuronNum))-Ypred_slc(slctr,timewindow(timepoint),NeuronNum);
                else
                pll_slc(tSilencingLag,NeuronNum) = pll_slc(tSilencingLag,NeuronNum) + ...
                                      mean(TrY(slctr,timewindow(timepoint),(NeuronNum))*log(Ypred_slc(:,timewindow(timepoint),NeuronNum))-Ypred_slc(:,timewindow(timepoint),NeuronNum));
                end
                  % this or above, and change slc timewindow

%                pll_ctrl(tSilencingLag,NeuronNum) = pll_ctrl(tSilencingLag,NeuronNum) + TrY(slctr,timewindow(timepoint),(NeuronNum))*log(mean(Ypred(:,timewindow(timepoint),NeuronNum),1))-mean(Ypred(:,timewindow(timepoint),NeuronNum),1);
%                pll_slc(tSilencingLag,NeuronNum) = pll_slc(tSilencingLag,NeuronNum) + TrY(slctr,timewindow(timepoint),(NeuronNum))*log(mean(Ypred_slc(:,timewindow(timepoint),NeuronNum),1))-mean(Ypred_slc(:,timewindow(timepoint),NeuronNum));
                pll_ctrl(tSilencingLag,NeuronNum) = pll_ctrl(tSilencingLag,NeuronNum) + ...
                    mean(TrY(slctr,timewindow(timepoint),(NeuronNum)).*log(Ypred(:,timewindow(timepoint),NeuronNum))-Ypred(:,timewindow(timepoint),NeuronNum));
             end
            end
%            mse_ctrl(tSilencingLag,NeuronNum) = mse_ctrl(tSilencingLag,NeuronNum)/size(TrY(:,timewindow,(NeuronNum)),1);
%            mse_slc(tSilencingLag,NeuronNum) = mse_slc(tSilencingLag,NeuronNum)/size(TrY(:,timewindow,(NeuronNum)),1);  
        end
        
    end
    
end
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
    for j = 1:8
        if j ~= i
            hold on; scatter(i,dtrue(i,:)*dhat(j,:)','c.')
            others(end+1) = dtrue(i,:)*dhat(j,:)';
            % hold on; scatter(i,dtrue(i,:)*dtrue(j,:)','g.')
        end
    end
end
xlim([0 9]);
subplot(1,2,2);scatter(randn(1,numel(self)),self,100,'k.');
hold on;scatter(10+randn(1,numel(others)),others,100,'c.');
xlim([-5,15])
% 1 and 4 don't work, I think needs noise measure
%figure;imagesc(dhat*dhat')
%figure;imagesc(dtrue*dtrue')
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
