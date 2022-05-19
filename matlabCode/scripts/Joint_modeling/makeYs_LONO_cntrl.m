function [Ypred,TrYcntrl] = makeYs_LONO_cntrl(res,UInd,fold)

Ncells = size(res.Xval.test_seq{fold}(1).y,1);

%%% these all from test trials of the particular fold
TrYcntrl = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,Ncells);
Ypred = nan(size(TrYcntrl));
try
    for neuron = 1:Ncells % held out neuron
        others = setdiff(1:Ncells,neuron); % All but the heldout neuron
        
        
        seqLM = []; % this is seq from all neurons but neuron, not necessarily LM        
        TrNumcntrl = []; % this will be identical for al neurons, only depends on the trials
        
        for tr = 1:numel(res.Xval.test_seq{fold})
            if sum(sum(res.Xval.test_seq{fold}(tr).u(2:end,:)))==0
                TrYcntrl(tr,:,neuron) = res.Xval.test_seq{fold}(tr).y(neuron,:)';
                TrNumcntrl(end+1) = tr;
                seqLM(end+1).y = res.Xval.test_seq{fold}(tr).y(others,:);
                seqLM(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
                seqLM(end).T = res.Xval.test_seq{fold}(tr).T;
            end
        end
        
        
        %%% get models from the traineing set, cut out irrelevant parts
        paramsLM = res.Xval.params{fold};
        paramsLM.model.C =  paramsLM.model.C(others,:);
        paramsLM.model.d=  paramsLM.model.d(others);
        
        % use params and test set activity of all but one neuron or
        % area to calculate the posterior
        [seqLM,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM);
        
        
        % use the parameters and posterior, to predict firing rate in
        % the left out neuron or area
        %     Ypred = nan(length(seqLM),res.Xval.train_seq{fold}(1).T,NV1Cells);
        %     Ypred_slc = nan(length(seqLM_slc),res.Xval.train_seq{fold}(1).T,NV1Cells);
        
        C = res.Xval.params{fold}.model.C(neuron,1:16);
        d = res.Xval.params{fold}.model.d(neuron);
        
        for tr = 1:length(seqLM)
            if 0
                Ypred{fold}(tr,:,:) = exp(C*seqLM(tr).posterior.xsm(1:16,:) + d)';
                %Ypred{fold}(tr,:,:) = (C*seqLM(tr).posterior.xsm(1:8,:) + d)';
            else
                Vsm   = reshape(seqLM(tr).posterior.Vsm' ,16,16,res.Xval.train_seq{fold}(1).T);
                for t = 1:res.Xval.train_seq{fold}(1).T
                    Sigma = Vsm(1:16,1:16,t); % are the first 2 correct?
                    Ypred(tr,t,neuron) = (exp(C*seqLM(tr).posterior.xsm(1:16,t) + d + diag(C*Sigma*C')/2)'); % check diag
                end
            end
            
        end
        
    end
    
    
    TrYcntrl = TrYcntrl(TrNumcntrl,:,:);
    Ypred = Ypred(find(~isnan(sum(sum(Ypred(:,:,:),2),3))),:,:);
catch
    TrNumcntrl = [];    
    for tr = 1:numel(res.Xval.test_seq{fold})
        if sum(sum(res.Xval.test_seq{fold}(tr).u(2:end,:)))==0
            TrYcntrl(tr,:,:) = res.Xval.test_seq{fold}(tr).y(:,:)';
            TrNumcntrl(end+1) = tr;
            
        end
    end
    TrYcntrl = TrYcntrl(TrNumcntrl,:,:);
    Ypred = nan(size(TrYcntrl));
    disp(['couldnt get states for this animal cntrl at fold ',num2str(fold)])
    
end
end