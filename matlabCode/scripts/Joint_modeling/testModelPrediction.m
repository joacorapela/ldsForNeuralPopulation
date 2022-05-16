animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

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

for animali  = 1%:6
    if CntrlOnly
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
    else
        cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
    end
    
  
    filelist = dir(['*splt1_onlyCorrect_exGrooming_go_5FoldXval*.mat']);
    res = load(filelist(end).name);
    nfold = 5;  
    
   Ypred = cell(1,nfold);
   Ypred_slc = cell(1,nfold);
   TrY = cell(1,nfold);
   TrYcntrl = cell(1,nfold);
    
    for    fold = 1:nfold
        
             
        tSilencingLag = 2;
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
        for tr = 1:numel(res.Xval.test_seq{fold})
            if sum(sum(res.Xval.test_seq{fold}(tr).u(2:end,:)))==0
                
                TrYcntrl{fold}(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(1:NV1Cells,:)';
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
                for t = 1:Xval.train_seq{fold}(1).T
                    Sigma = Vsm(1:8,1:8,t); % are the first 2 correct?
                    Ypred_slc{fold}(tr,t,:) = exp(C*seqLM_slc(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)';
                end
            end
        end
        
        for tr = 1:length(seqLM)
            if 0
                Ypred{fold}(tr,:,:) = exp(C*seqLM(tr).posterior.xsm(1:8,:) + d)';
            %Ypred{fold}(tr,:,:) = (C*seqLM(tr).posterior.xsm(1:8,:) + d)';
            else
                Vsm   = reshape(seqLM(tr).posterior.Vsm' ,16,16,res.Xval.train_seq{fold}(1).T);
                for t = 1:Xval.train_seq{fold}(1).T
                    Sigma = Vsm(1:8,1:8,t); % are the first 2 correct?
                    Ypred{fold}(tr,t,:) = exp(C*seqLM(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)'; % check diag
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
