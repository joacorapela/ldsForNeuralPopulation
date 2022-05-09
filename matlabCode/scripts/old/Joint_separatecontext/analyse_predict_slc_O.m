animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

orth = nan; % no orth here

allA = nan(6,16,16,2);
allB = nan(6,16,2);
for animali  = 1%:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
    
    %  d = dir(['*RND*_onlyCorrect_exGrooming_go.mat']);
    d = dir(['*RND7_onlyCorrect_exGrooming_go.mat']);
    
    rnd_vb = nan(1,length(d));
    for i = 1:length(d)
        if  numel(strfind(d(i).name,'onlyCorrect'))
            load(d(i).name);
            vb = varBound(~isnan(varBound));
            rnd_vb(i) = vb(end);
        end
    end
    [~,n]=max(rnd_vb);
    res = load(d(n).name);
       
     cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
     % lets say rnd
     slcres = load('Joint_PLDSfitRes_22_02_18_10_58_12_RND4_onlyCorrect_exGrooming_go.mat');
    
       NV1Cells = 27; % make dynamic
     % lag = 2;
     TrU = nan(numel(slcres.seq),slcres.seq(1).T);
     TrUv = nan(numel(slcres.seq),slcres.seq(1).T);
     TrY = nan(numel(slcres.seq),slcres.seq(1).T,NV1Cells);
     TrNum =[];
      TrYcntrl = nan(numel(slcres.seq),slcres.seq(1).T,NV1Cells);
     for tr = 1:numel(slcres.seq)
         if sum(slcres.seq(tr).u(2,:))==0%(slcres.seq(tr).u(2,62) - slcres.seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
         
         TrYcntrl(tr,:,:) =  slcres.seq(tr).y(1:NV1Cells,:)';
         elseif(slcres.seq(tr).u(2,62) - slcres.seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
        
             TrNum(end+1) = tr;
             TrU(tr,:) =  slcres.seq(tr).u(2,:);
             TrUv(tr,:) =  slcres.seq(tr).u(1,:);
             TrY(tr,:,:) =  slcres.seq(tr).y(1:NV1Cells,:)';
         end
     end
     TrU = TrU(TrNum,:);
     TrUv = TrUv(TrNum,:);
     TrY = TrY(TrNum,:,:);
     
     Slc_onset = 62;%find(nanmean(TrU(:,:)),1);
     nTimePoints = 12;
     
     nStates = 8; % 8
     
     % V1 states at the onset of silencing
    % V1X0 = arrayfun(@(x) x.posterior.xsm(1:nStates,Slc_onset),slcres.seq(TrNum),'UniformOutput',0); % cell length: numtrials
     V1X0 = arrayfun(@(x) x.posterior.xsm(1:nStates,Slc_onset),res.seq,'UniformOutput',0); % cell length: numtrials

     
     % 10 time points after onset:
     %    NV1Cells = 27; % make dynamic
     
     Ypred = nan(numel(TrNum),nTimePoints,NV1Cells);
     
     A = res.params.model.A(1:nStates,1:nStates); 
     B = res.params.model.B(1:nStates,:);
     C = res.params.model.C(1:NV1Cells,1:nStates);
     d = res.params.model.d(1:NV1Cells);
     for tr = 1:length(TrNum)
         X0 = V1X0{tr};
         
         for t = 1:nTimePoints
             Ypred(tr,t,:) = exp(C*X0 + d);
             X1 = A*X0+B*TrUv(tr,Slc_onset+t);% assumes u is 1 since during vis stimulus;
             X0 = X1;
         end
         
         
     end
     
     figure;
     NeuronList = [6 11 12 15 18 19 20 23 24 25];
     for NeuronNum = 1:10%1:27
         % trial averaged
         subplot(1,10,NeuronNum);plot(nanmean(Ypred(:,:,NeuronList(NeuronNum)),1))
         hold on;plot(nanmean(TrY(:,Slc_onset:Slc_onset+nTimePoints,NeuronList(NeuronNum)),1),'b')
         hold on;plot(nanmean(TrYcntrl(:,Slc_onset:Slc_onset+nTimePoints,NeuronList(NeuronNum)),1),'k')
         title(num2str(NeuronList(NeuronNum)))
         
     end
     
     
end
%% alternatively, give Y values and do inference with all states **** do res= slcres or not

animallist ={'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};

orth = nan; % no orth here

allA = nan(6,16,16,2);
allB = nan(6,16,2);
for animali  = 1%:6
    cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext_CntrlOnly/17msBins/'])
    
    %  d = dir(['*RND*_onlyCorrect_exGrooming_go.mat']);
    d = dir(['*RND0_onlyCorrect_exGrooming_go.mat']);
    
    rnd_vb = nan(1,length(d));
    for i = 1:length(d)
        if  numel(strfind(d(i).name,'onlyCorrect'))
            load(d(i).name);
            vb = varBound(~isnan(varBound));
            rnd_vb(i) = vb(end);
        end
    end
    [~,n]=max(rnd_vb);
    res = load(d(n).name);
       
     cd(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',animallist{animali},'/Joint_trial_based_splitContext/17msBins/'])
     % lets say rnd
     % slcres = load('Joint_PLDSfitRes_22_02_18_10_58_12_RND4_onlyCorrect_exGrooming_go.mat');
       d = dir(['*RND0_onlyCorrect_exGrooming_go.mat']);
        slcres = load(d(1).name);
    
       NV1Cells =27;% all 37 % make dynamic
       pvcell =[]%36
     % lag = 2;
     TrU = nan(numel(slcres.seq),slcres.seq(1).T);
     TrUv = nan(numel(slcres.seq),slcres.seq(1).T);
     TrY = nan(numel(slcres.seq),slcres.seq(1).T,NV1Cells);
     TrNum =[];
      TrYcntrl = nan(numel(slcres.seq),slcres.seq(1).T,NV1Cells);
      seqLM = [];

     for tr = 1:numel(slcres.seq)
         if sum(slcres.seq(tr).u(2,:))==0%(slcres.seq(tr).u(2,62) - slcres.seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
         
         TrYcntrl(tr,:,:) =  slcres.seq(tr).y(1:NV1Cells,:)';
        %  elseif(slcres.seq(tr).u(2,11) - slcres.seq(tr).u(2,10)) > 0 % change to else, if you want all silencing onsets not only lag 2
        % elseif(slcres.seq(tr).u(2,62) - slcres.seq(tr).u(2,61)) > 0 % change to else, if you want all silencing onsets not only lag 2
        elseif(slcres.seq(tr).u(2,85) - slcres.seq(tr).u(2,84)) > 0 
            seqLM(end+1).y =   slcres.seq(tr).y(NV1Cells+1:end,:);
            seqLM(end).y(pvcell-NV1Cells,:) = [];
             seqLM(end).u =   slcres.seq(tr).u(1,:); % 1 0r : if res= slcres
                  seqLM(end).T = slcres.seq(tr).T;
                  
             TrNum(end+1) = tr;
             TrU(tr,:) =  slcres.seq(tr).u(2,:);
             TrUv(tr,:) =  slcres.seq(tr).u(1,:);
             TrY(tr,:,:) =  slcres.seq(tr).y(1:NV1Cells,:)';
         end
     end
     TrU = TrU(TrNum,:);
     TrUv = TrUv(TrNum,:);
     TrY = TrY(TrNum,:,:);
     
     paramsLM = res.params;
     paramsLM.model.C(pvcell,:)=[];
     paramsLM.model.d(pvcell,:)=[];
     paramsLM.model.C =  paramsLM.model.C(NV1Cells+1:end,:);
     paramsLM.model.d=  paramsLM.model.d(NV1Cells+1:end);
     
     [SeqLM,~] = res.params.model.inferenceHandle(paramsLM,seqLM);
     
     Ypred = nan(numel(TrNum),res.seq(1).T,NV1Cells);
     Ypred_self = nan(numel(TrNum),res.seq(1).T,NV1Cells);
     
     C = res.params.model.C(1:NV1Cells,1:8);
     d = res.params.model.d(1:NV1Cells);
     for tr = 1:length(SeqLM)         
         Ypred(tr,:,:) = exp(C*SeqLM(tr).posterior.xsm(1:8,:) + d)';
         Ypred_self(tr,:,:) = exp(C*res.seq(tr).posterior.xsm(1:8,:) + d)';               
     end
     
     figure;
     NeuronList = [6 11 12 15 18 19 20 23 24 25];
     for NeuronNum = 1:10%1:27
         % trial averaged
         subplot(10,1,NeuronNum);plot(nanmean(Ypred(:,:,NeuronList(NeuronNum)),1),'c')
         hold on;plot(nanmean(Ypred_self(:,:,NeuronList(NeuronNum)),1),'r')
         hold on;plot(nanmean(TrY(:,:,NeuronList(NeuronNum)),1),'b')
         hold on;plot(nanmean(TrYcntrl(:,:,NeuronList(NeuronNum)),1),'k')
         title(num2str(NeuronList(NeuronNum)))
         xlim([50 100])
     end
     
     
end
