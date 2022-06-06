% cd (rootdir)
% for animali = animalIds
%     % load the lono/loao file
%     filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
%     res = load(fullfile(rootdir,filename));
%     
%     if CntrlOnly
%         cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
%         filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
%     else
%         cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
%         filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
%     end
%     mdl = load(filelist(end).name);
%     TrYcntrl = res.TrYcntrl;% not useful, only v1
%     vecs = [];
%     cntrltrl = [];
%     fold = 3;
%     for tr = 1:numel(mdl.Xval.train_seq{fold})
%         if sum(sum(mdl.Xval.train_seq{fold}(tr).u(2:end,:)))==0
%             cntrltrl(end+1) = tr;
%         end
%     end
%     for tSilencingLag = 1:8
%         
%         Ypred_slc = res.Lag{tSilencingLag}.Ypred_slc;
%         TrY = res.Lag{tSilencingLag}.TrY;
%         
%         if ~CntrlOnly
%             timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
%         else
%             timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
%         end
%         timewindow = timewindow(1:end);
%         
%         %vecs = [vecs; squeeze(nanmean(TrYcntrl(:,timewindow,:),2))]; % trials by neurons
%         % wrong though, we want LM's activity, not V1
%         % vecs = [vecs; nanmean(squeeze(nanmean(TrYcntrl(:,timewindow,:),2)),1)]; % trials by neurons
%         LMneuronsByTrials = cell2mat(arrayfun(@(x) nanmean(x.y(res.NV1Cells+1:end,timewindow),2),mdl.Xval.train_seq{fold},'UniformOutput',0));
%         vecs = [vecs;nanmean(LMneuronsByTrials(:,cntrltrl),2)']; % this should be LM, from one fold, also, onlyfromcntrl trials
%         %vecs = [vecs;LMneuronsByTrials(:,cntrltrl)'];
%     end
%     figure
%     for i=1:4
%         for fold = 1:5
%             A = mdl.Xval.params{fold}.model.A;
%             [CO,W] = orth_c_svd(mdl.Xval.params{fold});
%             A = W*A*(W');
%             [U,S,V] = svds(A(1:8,9:16));
%             [c,s,l] = pca(vecs,'Centered',0); % zero or one?
%             
%             sv = V(:,i); % V or L?
%             %projsvLM = (mdl.Xval.params{fold}.model.C(res.NV1Cells+1:end,9:16)*sv+mdl.Xval.params{fold}.model.d(res.NV1Cells+1:end))';
%             projsvLM = (CO(res.NV1Cells+1:end,9:16)*sv+mdl.Xval.params{fold}.model.d(res.NV1Cells+1:end))';
%             projsvLM = projsvLM/norm(projsvLM);
%             subplot(1,4,i);hold on;plot(abs(projsvLM*c))
%         end
%     end
% end
% 
% %%
% cd (rootdir)
% CntrlOnly = 0;
% 
% %%% orth necessary? PCA with same bin size? Av or pool? pre subtract or
% %%% not? Train or test per fold?
% %%%
% figure
% for animali = animalIds
%     % load the lono/loao file
%     filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
%     res = load(fullfile(rootdir,filename));
%     
%     if CntrlOnly
%         cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
%         filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
%     else
%         cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
%         filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
%     end
%     mdl = load(filelist(end).name);
% 
%     Win = 58:116*3/4;%58:100;%58:116*3/4;%58:116*3/4; % viaul window
%     pre = 1:57;
% 
%     for fold = 1:5 % for now, it is train, not test
%         cntrltrl = [];
%         for tr = 1:numel(mdl.Xval.train_seq{fold})
%             if sum(sum(mdl.Xval.train_seq{fold}(tr).u(2:end,:))) == 0
%                 cntrltrl(end+1) = tr;
%             end
%         end
%         LMAcStates = arrayfun(@(x) x.posterior.xsm(9:16,:),mdl.Xval.train_seq{fold}(cntrltrl),'UniformOutput',0); % cell of nTr elements, each 8 *116
%        
%         %Av = zeros(8,numel(Win));
%         Av = [];
%         for tr = 1:numel(LMAcStates)
%             Av = cat(3,Av,LMAcStates{tr}(:,Win) - nanmean(LMAcStates{tr}(:,pre),2));
%            % Av = Av +LMAcStates{tr}(:,Win);
%            %  Av = Av +(LMAcStates{tr}(:,Win) - nanmean(LMAcStates{tr}(:,pre),2));
%         end
%        % Av = Av/numel(LMAcStates);
%        % Av = nanmean(Av,3); 
%        Av = reshape(Av,8,[]);
%        
%         A = mdl.Xval.params{fold}.model.A;
%         %            [CO,W] = orth_c_svd(mdl.Xval.params{fold});
%         %            A = W*A*(W');
%         [U,S,V] = svds(A(1:8,9:16));
%         [c,s,l] = pca(Av','Centered',0); % zero or one?
%         for i=1:4
%             sv = V(:,i)'; % V or L?            
%             subplot(6,4,4*(animali-1)+i);hold on;plot(abs(sv*c))
%         end
%     end
% end

% cd (rootdir)
% CntrlOnly = 0;
% for avOrPool = 1:2
%     for subPre = 0:1
%         
%         for useSlcWin = 1% 0:1
%             for doCenter = 0%0:1
%                 for winOpt = 1% doesn't matter:3
%                     for acOrSt = 2%1:2 
%                         for testOrTrain = 1% reain more reliable version1:2
%                             try
%                                 draft_svAlignment
%                             catch
%                                 disp('caught')
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
%% move to test trials
%%% orth necessary? PCA with same bin size? Av or pool? pre subtract or
%%% not? Train or test per fold? Win or singe bins?
%%%
%%% TODO: add average per panel
% best: 2 - 0
avOrPool = 2; % 1 av 2 pool
subPre = 0 ; 

% Fix:
useSlcWin = 1; % if 1: tw if 0: Win
doCenter = 0;
winOpt = 1; % 1,2,3 - doesnt matter
acOrSt = 2;%2 state ,1 activity
testOrTrain = 1;% test 1 train 2
%maybe add longrange or local

runID = [num2str(avOrPool),num2str(subPre),num2str(useSlcWin),num2str(doCenter),num2str(winOpt),num2str(acOrSt),num2str(testOrTrain)];
    
pre = 1:57; 
if winOpt == 1
    Win = 58:116*3/4;
elseif winOpt == 2
    Win = 1:116;
elseif winOpt == 3
    Win = 58:100;
end
f = figure;
for animali = animalIds
    % load the lono/loao file
    filename = ['animal',num2str(animali),'_',leaveout,'_','CntrlOnly',num2str(CntrlOnly)];
    res = load(fullfile(rootdir,filename));
    tw = cell(1,8);
    for tSilencingLag = 1:8   
        if ~CntrlOnly
            timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.train_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
        else
            timewindow = find(sum(cell2mat(cellfun(@(y) y(tSilencingLag+1,:), arrayfun(@(x) x.u, mdl.Xval.test_seq{1},'UniformOutput',0),'UniformOutput',0)'),1)>1);
        end
        tw{tSilencingLag} = timewindow(1:end);
    end
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
        filelist = dir(['*splt1_RND0_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
        filelist = dir(['*splt1_onlyCorrect_exGrooming_go_',num2str(nfold),'FoldXval*.mat']);
    end
    mdl = load(filelist(end).name);
    
    for fold = 1:5 %
        if testOrTrain == 1
            [mdl.Xval.test_seq{fold},~] = mdl.Xval.params{fold}.model.inferenceHandle(mdl.Xval.params{fold},mdl.Xval.test_seq{fold});
            cntrltrl = [];
            for tr = 1:numel(mdl.Xval.test_seq{fold})
                if sum(sum(mdl.Xval.test_seq{fold}(tr).u(2:end,:))) == 0
                    cntrltrl(end+1) = tr;
                end
            end
            if acOrSt == 2
                LMAcStates = arrayfun(@(x) x.posterior.xsm(9:16,:),mdl.Xval.test_seq{fold}(cntrltrl),'UniformOutput',0); % cell of nTr elements, each 8 *116
            elseif acOrSt == 1
                LMAcStates = arrayfun(@(x) x.y(res.NV1Cells+1:end,:),mdl.Xval.test_seq{fold}(cntrltrl),'UniformOutput',0);
            end
        elseif testOrTrain == 2
            cntrltrl = [];
            for tr = 1:numel(mdl.Xval.train_seq{fold})
                if sum(sum(mdl.Xval.train_seq{fold}(tr).u(2:end,:))) == 0
                    cntrltrl(end+1) = tr;
                end
            end
            if acOrSt == 2
                LMAcStates = arrayfun(@(x) x.posterior.xsm(9:16,:),mdl.Xval.train_seq{fold}(cntrltrl),'UniformOutput',0); % cell of nTr elements, each 8 *116
            elseif acOrSt == 1
                LMAcStates = arrayfun(@(x) x.y(res.NV1Cells+1:end,:),mdl.Xval.train_seq{fold}(cntrltrl),'UniformOutput',0);
            end
        end
        
        
        Av = [];
        for tr = 1:numel(LMAcStates)
            if useSlcWin
                Av = cat(3,Av,cell2mat(cellfun(@(x) nanmean(LMAcStates{tr}(:,x),2),tw,'UniformOutput',0)));
            else
                if subPre
                    Av = cat(3,Av,LMAcStates{tr}(:,Win) - nanmean(LMAcStates{tr}(:,pre),2));
                else
                    Av = cat(3,Av,LMAcStates{tr}(:,Win));
                end
            end
        end
        if avOrPool == 1
            Av = nanmean(Av,3);
        elseif avOrPool == 2
            Av = reshape(Av,8,[]);
        end

        A = mdl.Xval.params{fold}.model.A;
        Alv = A(1:8,9:16);
        %            [CO,W] = orth_c_svd(mdl.Xval.params{fold});
        %            A = W*A*(W');
        [U,S,V] = svd(Alv);
        [c,s,l] = pca(Av','Centered',doCenter); % zero or one?
        for i=1:4           
            if acOrSt == 2
                sv = V(:,i)';
                subplot(6,4,4*(animali-1)+i);hold on;plot(abs(sv*c),'color',[.6 .6 .6]);
            elseif acOrSt == 1
                sv = V(:,i)';
                projsvLM = (mdl.Xval.params{fold}.model.C(res.NV1Cells+1:end,9:16)*sv'+mdl.Xval.params{fold}.model.d(res.NV1Cells+1:end))';
                projsvLM = projsvLM/norm(projsvLM);
                subplot(6,4,4*(animali-1)+i);hold on;plot(abs(projsvLM*c))                 
            end
            
        end
    end
end
f.Name = runID;
