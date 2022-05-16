function [seq,N_V1,N_LM] = buildTrialBasedSeq_JointSepCon(summarymatfile, binSizems,binWinms,splitDelays,perffile,trialType,Xval,CntrlOnly)
res = load(fullfile(summarymatfile.folder,summarymatfile.name));
if ~isnan(binWinms)
    error('custized window for trial, not implemented yet')
end

% not implemented yet for joint:
% if LONO.do
%     LONO.fileContents = load(fullfile(LONO.file.folder,LONO.file.name));
%     res = assignToFold(LONO.fileContents.Fold{LONO.fold}.trainInd,res);
% end

if 1 % condition on trialtype, keep lono 0
    if ~isempty(perffile)
        p = load(fullfile(summarymatfile.folder,perffile.name));
        res.correctgotrialind = p.correctgotrialind;
        res.correctnogotrialind = p.correctnogotrialind;
    end
    
    res.correcttrialind = [res.correctgotrialind,res.correctnogotrialind];
    
    if strcmp(trialType,'onlyCorrect_exGrooming')
        validTrIdx = intersect(res.correcttrialind,res.nogroomingind);
    elseif strcmp(trialType,'onlyCorrect_exGrooming_go')
        validTrIdx = intersect(res.correctgotrialind,res.nogroomingind);
    elseif strcmp(trialType,'onlyCorrect_exGrooming_nogo')
        validTrIdx = intersect(res.correctnogotrialind,res.nogroomingind);
    elseif strcmp(trialType,'exGrooming')
        validTrIdx = res.nogroomingind;
    else
        error('trial type not implemented')
    end
    res = assignToFold(validTrIdx,res);
end

units_V1 = [res.V1{1}.SingleUnitSpikeTimes,res.V1{2}.SingleUnitSpikeTimes]; % cell of length Nneurons,
N_V1 = length(units_V1); % number of neurons                                                                                                    % containing spike times

units_LM = [res.LM{1}.SingleUnitSpikeTimes,res.LM{2}.SingleUnitSpikeTimes]; % cell of length Nneurons,
N_LM = length(units_LM);

units = [units_V1,units_LM];

LaserTrials = find(~isnan(res.LaserDelay));

edgestep = binSizems*30; % binsize in samples
nTr = size(res.PAllOn,1); % the number of all trials

LaserDur = 150*30; % length of laser stimulation in samples, add shape later. For now, fixed.

% pre allocating the struct array: makes it much faster
for Tr = size(res.PAllOn,1):-1:1
    seq(Tr) = struct;
end


if  splitDelays
    % error('splitDelay 0 not implemented')
    nInputs = 9;
else
    nInputs =2;
end
res.allLaserD = makeDiscreteDelayInd(res.LaserDelayBinned);

Discard = [];
noLaser = [];
for Tr = 1:size(res.PAllOn,1)
    % res.PAllOn(Tr,:) contains the time stamps of the samples in trial Tr from -1sec to 1sec
    
    middlep = ceil(size(res.PAllOn,2)/2);
    edges = [fliplr(res.PAllOn(Tr,middlep-edgestep):-edgestep:res.PAllOn(Tr,1)),res.PAllOn(Tr,middlep):edgestep:res.PAllOn(Tr,end)];
    T = size(edges,2) - 1; % number of bins
    
    % spikes
    seq(Tr).y = nan((N_V1+N_LM),T);
    for neuron = 1:(N_V1+N_LM)
        [seq(Tr).y(neuron,:),~] = histcounts(units{neuron},edges);
    end
    
    % make u matrix: order: vg,vng,lg,lng
    seq(Tr).u = zeros(nInputs,T);
    binStart = edges(1:T);
    binEnd = edges(2:end);
    
    % visual stimulus
    pdtrace = zeros(1,T); % photo diode signal
    pdtrace(find(binEnd>res.PStepTimeStampOn(Tr) & binEnd<res.PStepTimeStampOff(Tr))) = 1;
    
    % laser stimulus
    lstrace = zeros(1,T);
    LaserTrInd = find(LaserTrials == Tr);
    if numel(LaserTrInd) % if trial has laser stimulation
        lstrace(find(binEnd>res.LStepTimeStampOn(LaserTrInd) & binEnd<(res.LStepTimeStampOn(LaserTrInd)+LaserDur))) = 1;
    end
    
    % find laserD, if not splitting input, =0, otherwise value in [0-7]
    % depending on laser delay wrt the stimulus onset
    % doublecheck that nans are aligned with no laser
    if ~splitDelays
        laserD = 0;
    else
        laserD = res.allLaserD(Tr); % could be nan (no laser), 0-7 Delays, or 8 : dicard the trial
        if isnan(laserD)
            laserD = 0;
            noLaser = [noLaser,Tr];
        elseif laserD == 8
            lstrace = nan(size(lstrace));
            laserD = 0;
            Discard = [Discard,Tr];
        end
        % TODO: assign delay values in ms here as well. (in a separate
        % array)
    end
    
    if 1%~numel(LaserTrInd) % other option for when only using cntrl trials
        seq(Tr).u(1,:) = pdtrace;
        seq(Tr).u(2+laserD,:) = lstrace;
    else
        Discard = [Discard,Tr];
    end
    
    seq(Tr).T = T;
end

if splitDelays
    if any(find(arrayfun(@(x) sum(sum(x.u(2:end,:))) == 0,seq,'UniformOutput',1)) == noLaser) == 0
        error('error in extracting laser delays. No laser trials dontmatch up')
    end
end

seq(Discard) = [];
res.allLaserD(Discard) = [];
%%%% divide into training and test sets
if ~CntrlOnly
    
    if Xval.do
        
        Xval.train_seq = cell(1,Xval.nFolds);
        Xval.test_seq = cell(1,Xval.nFolds);
        
        if strcmp(Xval.method,'rand')
            rng(Xval.Seed,'twister')
            indices = crossvalind('Kfold',numel(seq),Xval.nFolds);
        elseif strcmp(Xval.method,'lag') & ~splitDelays
            indices = res.allLaserD; % contains either 0-7 for lags or nan for control trials
            % also, divide nan values (control trials) in between the
            % folds, randomly [test contains 1 lag + 1/8 of control trials]
            CnrlIndices = find(isnan(indices));
            rng(Xval.Seed,'twister')
            indices(CnrlIndices) = randi(8,[1,numel(CnrlIndices)])-1;
            % bringing the 0-7 range to 1-8 for fold assignment. The
            % unassigned lasers are not discarded when ~splitDelays. But
            % they will get a value of 9, and will never be in test, always
            % in training set
            indices = indices + 1;
        elseif strcmp(Xval.method,'lag') & splitDelays
            error('lag method cant work when splitting delays (test and train models different)')    
        end
        
        for FoldN = 1:Xval.nFolds
            Xval.Fold{FoldN}.testInd = (find(indices == FoldN))';
            Xval.Fold{FoldN}.trainInd = (find(~(indices == FoldN)))';
            
            Xval.train_seq{FoldN} = seq(Xval.Fold{FoldN}.trainInd);
            Xval.test_seq{FoldN} = seq(Xval.Fold{FoldN}.testInd);
        end
        
        % Xval.train_seq{1} to access previous seq
        seq = Xval;
    end
else
    if Xval.do
        if  strcmp(Xval.method,'lag')
            error('lag method doesnt exist for CntrlOnly, change to rand')
        end
        
        seqcntrl = seq;
        rest = seq;
        for tr = numel(seq):-1:1
            if sum(sum(seq(tr).u(2:end,:))) ~= 0
                seqcntrl(tr)=[];
            else
                rest(tr)=[];
            end
        end
        
        
        Xval.train_seq = cell(1,Xval.nFolds);
        Xval.test_seq = cell(1,Xval.nFolds);
        
        rng(Xval.Seed,'twister')
        indices = crossvalind('Kfold',numel(seqcntrl),Xval.nFolds);
        
        for FoldN = 1:Xval.nFolds
            Xval.Fold{FoldN}.testInd = (find(indices == FoldN))';
            Xval.Fold{FoldN}.trainInd = (find(~(indices == FoldN)))';
            
            Xval.train_seq{FoldN} = seqcntrl(Xval.Fold{FoldN}.trainInd);
            Xval.test_seq{FoldN} = seqcntrl(Xval.Fold{FoldN}.testInd);
            
            Xval.test_seq{FoldN} = [Xval.test_seq{FoldN},rest];
            
            % clear u column of control trials
            for tr = 1:numel(Xval.train_seq{FoldN})
                Xval.train_seq{FoldN}(tr).u(2:end,:) = [];
            end
            
            
        end
        
        
        seq = Xval;
    else
        
        for tr = numel(seq):-1:1
            if sum(sum(seq(tr).u(2:end,:))) ~= 0
                seq(tr)=[];
            else
                seq(tr).u(2:end,:) = [];
            end
        end
        
    end
end

% double check: ind(arrayfun(@(x)
% numel(find(isnan(x.u))),seq,'UniformOutput',1)) = Discard