
function trial_ll = heldout_loglike(resTr,LONOparams,FittedFold,heldoutN,doPlot)
% make y of the test dataset, leaving out heldoutN
% seq_one will only get a y_orig field which is its real spiking rate
[seq_minusOne,seq_One] = buildTrialBasedSeq_heldout(resTr.config.summarymatfile, resTr.config.binSizems,resTr.config.binWinms,...
    resTr.config.area,LONOparams.Fold{FittedFold}.testInd,heldoutN,resTr.config.splitDelays);

% inference: get x - before here, prepare the loading matrix (C) and d : eliminate the corresponding row)
% check if anything else is needed
resTs_minusOne = resTr;
resTs_minusOne.params.model.C(heldoutN,:) = [];
resTs_minusOne.params.model.d(heldoutN) = [];
resTs_minusOne.seq = seq_minusOne;
[resTs_minusOne.seq,~] = resTr.params.model.inferenceHandle(resTs_minusOne.params,resTs_minusOne.seq);
%                                                                                                                                                                  % (it has now a field named yorig with one row, no u)

pred = nan(length(seq_One),resTs_minusOne.seq(1).T); % trials*timebins
yOrig = nan(length(seq_One),resTs_minusOne.seq(1).T); % trials*timebins
for tr = 1:length(seq_One)
    % for trial tr:
    z = resTr.params.model.C(heldoutN,:) * resTs_minusOne.seq(tr).posterior.xsm + ...
        resTr.params.model.d(heldoutN);
    pred(tr,:) = exp(z);
    yOrig(tr,:) = seq_One(tr).yOrig;
end
% pred is actually the expectation of the poisson distribution
% gamma is hard coded here
if doPlot
    figure;plot(nanmean(pred,1));hold on;plot(nanmean(yOrig,1));legend({'prediction','original trace'})
end
% compare pred and yOrig, poisson loglikelihood: (double check)
trial_ll = nan(1,length(seq_One));
for tr = 1:length(seq_One)
    % for trial tr:
    lambda = pred(tr,:);
    x = yOrig(tr,:);
    n = length(yOrig(tr,:));
    
    % this is a 1*n array, each element loglikelihood of one timepoint
    ll = -lambda + x.*log(lambda) - log(factorial(x));
    % ll of individual time bins are summed to make the total ll of the
    % trial
    trial_ll(tr) = sum(ll);
    
end
end

