function seq = AddinputDyn(seq)


stimInd = find(sum(abs(seq(1).u(1:2,:)),1));
stimLength = length(stimInd);
trialLength = seq(1).T;

for T=1:length(seq)
    fullU = seq(T).u;
    non_visual_u = fullU(3:end,:);
    
    goStack = zeros(stimLength,trialLength);
    nogoStack = zeros(stimLength,trialLength);
    
    
    for i = 1:stimLength
        goStack(i,stimInd(i)) = fullU(1,stimInd(i));
        nogoStack(i,stimInd(i)) = fullU(2,stimInd(i));
    end
    
    seq(T).u = [goStack;nogoStack;non_visual_u];
    
end