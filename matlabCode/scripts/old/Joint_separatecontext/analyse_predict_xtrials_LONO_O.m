animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

orth = nan; % no orth here
CntrlOnly = 1; % 0 or not


if CntrlOnly
    UInd = 1;
else
    UInd = 1:2;
end

for animali  = 1%:6
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    end
    
    %  d = dir(['*RND*_onlyCorrect_exGrooming_go.mat']);
    filelist = dir(['*RND0_onlyCorrect_exGrooming_go_Xval*.mat']);
    filelist = filelist(end-4:end);
    nfold = 5;
    
    
    Ypred = cell(1,nfold);
    Ypred_slc = cell(1,nfold);
    TrY = cell(1,nfold);
    TrYcntrl = cell(1,nfold);
    
    for    fold = 1:nfold
        
        res = load(filelist(fold).name);
        
        
        NV1Cells =11;%find(sum(res.Xval.params{1}.model.C(:,1:8),2)==0,1)-1%27;% all 37 % make dynamic maybe try animals with more lm cells
        Ncells = size(res.Xval.params{1}.model.C,1);
        restcells = setdiff(1:Ncells,NV1Cells);
        pvcell =[]%36
        % lag = 2;
        
        %%% these all from test trials
        TrU = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
        TrUv = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
        TrY{fold} = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,1);
        TrYcntrl{fold} = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,1);
        seqLM_slc = [];  % from test trials
        seqLM = [];
        TrNum =[];
        for tr = 1:numel(res.Xval.test_seq{fold})
            if sum(res.Xval.test_seq{fold}(tr).u(2,:))==0%(slcres.seq(tr).u(2,62) - slcres.seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
                
                TrYcntrl{fold}(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(NV1Cells,:)';
                seqLM(end+1).y = res.Xval.test_seq{fold}(tr).y(restcells,:);
                %seqLM(end).y(pvcell-NV1Cells,:) = [];
                seqLM(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
                seqLM(end).T = res.Xval.test_seq{fold}(tr).T;
                %  elseif(slcres.seq(tr).u(2,11) - slcres.seq(tr).u(2,10)) > 0 % change to else, if you want all silencing onsets not only lag 2
            elseif(res.Xval.test_seq{fold}(tr).u(2,62) - res.Xval.test_seq{fold}(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
            %elseif(res.Xval.test_seq{fold}(tr).u(2,85) - res.Xval.test_seq{fold}(tr).u(2,84)) > 0
                seqLM_slc(end+1).y = res.Xval.test_seq{fold}(tr).y(restcells,:);
                %seqLM_slc(end).y(pvcell-NV1Cells,:) = [];
                seqLM_slc(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
                seqLM_slc(end).T = res.Xval.test_seq{fold}(tr).T;
                
                TrNum(end+1) = tr;
                TrU(tr,:) =  res.Xval.test_seq{fold}(tr).u(2,:);
                TrUv(tr,:) =  res.Xval.test_seq{fold}(tr).u(1,:);
                TrY{fold}(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(NV1Cells,:)';
            end
        end
        TrU = TrU(TrNum,:);
        TrUv = TrUv(TrNum,:);
        TrY{fold} = TrY{fold}(TrNum,:,:);
        
        paramsLM = res.Xval.params{fold};
        paramsLM.model.C(pvcell,:)=[];
        paramsLM.model.d(pvcell,:)=[];
        paramsLM.model.C =  paramsLM.model.C(restcells,:);
        paramsLM.model.d=  paramsLM.model.d(restcells);
        
        [seqLM_slc,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM_slc);
        [seqLM,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM);
        
        Ypred{fold} = nan(length(seqLM),res.Xval.train_seq{fold}(1).T,1);
        Ypred_slc{fold} = nan(length(seqLM_slc),res.Xval.train_seq{fold}(1).T,1);
        
        C = res.Xval.params{fold}.model.C(NV1Cells,1:8);
        d = res.Xval.params{fold}.model.d(NV1Cells);
        for tr = 1:length(seqLM_slc)
            Ypred_slc{fold}(tr,:,:) = exp(C*seqLM_slc(tr).posterior.xsm(1:8,:) + d)';
        end
        
        for tr = 1:length(seqLM)
            Ypred{fold}(tr,:,:) = exp(C*seqLM(tr).posterior.xsm(1:8,:) + d)';
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

    
    
    figure;
    NeuronList =1:10;% [5 11 12 15 18 19 20 23 24 25];% 34 for animal 6, 22 qnd 24 for animal 3
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
