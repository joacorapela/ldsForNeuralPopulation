function  out = assignToFold(tInd,res)

% assign values
out.LaserDelay = res.LaserDelay(tInd);
out.LaserDelayBinned = res.LaserDelayBinned(tInd);
out.PAllOn = res.PAllOn(tInd,:);
out.PStepTimeStampOn = res.PStepTimeStampOn(tInd);
out.PStepTimeStampOff = res.PStepTimeStampOff(tInd);
out.licks = res.licks(tInd,:);
% test and double check this:
LaserTrialInd = find(~isnan(res.LaserDelay));
[~,TLaserTrialInd] = intersect(LaserTrialInd,tInd);
out.LStepTimeStampOn = res.LStepTimeStampOn(TLaserTrialInd);
out.LAllOn = res.LAllOn(TLaserTrialInd,:);
% and no laser off field(ramp down)
%
[~,~,out.gotrialind]= intersect(res.gotrialind,tInd);
[~,~,out.nogotrialind] = intersect(res.nogotrialind,tInd);
[~,~,out.correctgotrialind] = intersect(res.correctgotrialind,tInd);
[~,~,out.correctnogotrialind] = intersect(res.correctnogotrialind,tInd);
[~,~,out.nogroomingind] = intersect(res.nogroomingind,tInd);

%%% pass on timeseries defined continuously, not in trial
%%% structure
out.V1 = res.V1;
out.LM = res.LM;

