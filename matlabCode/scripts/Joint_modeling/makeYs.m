function [Ypred,Ypred_slc,TrYcntrl,TrY] = makeYs(res,NV1Cells,pvcell,UInd,fold,tSilencingLag)

%%% these all from test trials of the particular fold
TrU = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
TrUv = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T);
TrY = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,NV1Cells);
TrYcntrl = nan(numel(res.Xval.test_seq{fold}),res.Xval.test_seq{fold}(1).T,NV1Cells);
seqLM_slc = [];  % from test trials
seqLM = [];
TrNum =[];
TrNumcntrl = [];
for tr = 1:numel(res.Xval.test_seq{fold})
    if sum(sum(res.Xval.test_seq{fold}(tr).u(2:end,:)))==0
        TrYcntrl(tr,:,:) = res.Xval.test_seq{fold}(tr).y(1:NV1Cells,:)';
        TrNumcntrl(end+1) = tr;
        seqLM(end+1).y = res.Xval.test_seq{fold}(tr).y(NV1Cells+1:end,:);
        seqLM(end).y(pvcell-NV1Cells,:) = [];
        seqLM(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
        seqLM(end).T = res.Xval.test_seq{fold}(tr).T;
        
    elseif sum(sum(res.Xval.test_seq{fold}(tr).u(1+tSilencingLag,:))) > 0
        seqLM_slc(end+1).y = res.Xval.test_seq{fold}(tr).y(NV1Cells+1:end,:);
        seqLM_slc(end).y(pvcell-NV1Cells,:) = [];
        % in CntrlOnly, since there is no nonvisual B, non visual us of
        % test set is also removed
        seqLM_slc(end).u = res.Xval.test_seq{fold}(tr).u(UInd,:); % 1 0r : if res= slcres
        seqLM_slc(end).T = res.Xval.test_seq{fold}(tr).T;
        
        TrNum(end+1) = tr;
        TrU(tr,:) =  res.Xval.test_seq{fold}(tr).u(1+tSilencingLag,:); %% check this, modified
        %    TrU(tr,:) =  res.Xval.test_seq{fold}(tr).u(2,:); %% check this, modified
        TrUv(tr,:) =  res.Xval.test_seq{fold}(tr).u(1,:);
        TrY(tr,:,:) =  res.Xval.test_seq{fold}(tr).y(1:NV1Cells,:)';
    end
end


TrU = TrU(TrNum,:);
TrUv = TrUv(TrNum,:);
TrY = TrY(TrNum,:,:);
TrYcntrl = TrYcntrl(TrNumcntrl,:,:);

%%% get models from the traineing set, cut out irrelevant parts
paramsLM = res.Xval.params{fold};
paramsLM.model.C(pvcell,:)=[];
paramsLM.model.d(pvcell,:)=[];
paramsLM.model.C =  paramsLM.model.C(NV1Cells+1:end,:);
paramsLM.model.d=  paramsLM.model.d(NV1Cells+1:end);

% use params and test set activity of all but one neuron or
% area to calculate the posterior
[seqLM_slc,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM_slc);
[seqLM,~] = res.Xval.params{fold}.model.inferenceHandle(paramsLM,seqLM);


% use the parameters and posterior, to predict firing rate in
% the left out neuron or area
Ypred = nan(length(seqLM),res.Xval.train_seq{fold}(1).T,NV1Cells);
Ypred_slc = nan(length(seqLM_slc),res.Xval.train_seq{fold}(1).T,NV1Cells);

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
            Ypred_slc(tr,t,:) = (exp(C*seqLM_slc(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)');
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
            Ypred(tr,t,:) = (exp(C*seqLM(tr).posterior.xsm(1:8,t) + d + diag(C*Sigma*C')/2)'); % check diag
        end
    end
    
end


end