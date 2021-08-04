function [seq] = buildTrialBasedSeq(summarymatfile, binSizems,binWinms,area,LONO)
res = load(fullfile(summarymatfile.folder,summarymatfile.name));
if ~isnan(binWinms)
    error('custized window for trial, not implemented yet')
end

if LONO.do 
    res = res.Fold{LONO.fold}.train;
end


units = eval(sprintf('[res.%s{1}.SingleUnitSpikeTimes,res.%s{2}.SingleUnitSpikeTimes]',area,area)); % cell of length Nneurons, 
N = length(units); % number of neurons                                                                                                    % containing spike times
LaserTrials = find(~isnan(res.LaserDelay));

edgestep = binSizems*30; % binsize in samples
nTr = size(res.PAllOn,1); % the number of all trials

LaserDur = 150*30; % length of laser stimulation in samples, add shape later. For now, fixed.

% pre allocating the struct array: makes it much faster
for Tr = size(res.PAllOn,1):-1:1
    seq(Tr) = struct;
end

for Tr = 1:size(res.PAllOn,1)
    % res.PAllOn(Tr,:) contains the time stamps of the samples in trial Tr from -1sec to 1sec
    
    middlep = ceil(size(res.PAllOn,2)/2);                      
    edges = [fliplr(res.PAllOn(Tr,middlep-edgestep):-edgestep:res.PAllOn(Tr,1)),res.PAllOn(Tr,middlep):edgestep:res.PAllOn(Tr,end)];
    T = length(edges) - 1; % number of bins
    
    % spikes
    seq(Tr).y = nan(N,T);
    for neuron = 1:N
        [seq(Tr).y(neuron,:),~] = histcounts(units{neuron},edges);
    end
    
    % make u matrix: order: vg,vng,lg,lng
    seq(Tr).u = zeros(4,T); 
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
    
     % assign pdtrace to go or nogo trace
    if numel(find(res.gotrialind==Tr))
        seq(Tr).u(1,:) = pdtrace;
        seq(Tr).u(3,:) = lstrace;
    else
        seq(Tr).u(2,:) = pdtrace;
        seq(Tr).u(4,:) = lstrace;
    end
    
    seq(Tr).T = T;
end
