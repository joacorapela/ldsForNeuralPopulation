
function trial_ll = test_loglike(resTr,LONOparams,FittedFold,doPlot)

% make y of the test dataset
seq = buildTrialBasedSeq_test(resTr.config.summarymatfile, resTr.config.binSizems,resTr.config.binWinms,...
    resTr.config.area,LONOparams.Fold{FittedFold}.testInd,resTr.config.splitDelays);

% inference: get x
resTs = resTr;
resTs.seq = seq;
[resTs.seq,~] = resTr.params.model.inferenceHandle(resTs.params,resTs.seq);
%                                                                                                                                                                  % (it has now a field named yorig with one row, no u)
num_neurons = size(resTs.seq(1).y,1);

trial_ll = nan(num_neurons,length(seq)); % num neurons * num trials
for Neuron = 1:num_neurons
    
    pred = nan(length(seq),resTs.seq(1).T); % trials*timebins
    yOrig = nan(length(seq),resTs.seq(1).T); % trials*timebins
    for tr = 1:length(seq)
        % for trial tr:
        z = resTr.params.model.C(Neuron,:) * resTs.seq(tr).posterior.xsm + ...
            resTr.params.model.d(Neuron);
        pred(tr,:) = exp(z);
        yOrig(tr,:) = seq(tr).y(Neuron,:);
    end
    % pred is actually the expectation of the poisson distribution
    % gamma is hard coded here
    
    % compare pred and yOrig, poisson loglikelihood: (double check)
    
    for tr = 1:length(seq)
        % for trial tr:
        lambda = pred(tr,:);
        x = yOrig(tr,:);
        n = length(yOrig(tr,:));
        
        % this is a 1*n array, each element loglikelihood of one timepoint
        ll = -lambda + x.*log(lambda) - log(factorial(x));
        % ll of individual time bins are summed to make the total ll of the
        % trial
        trial_ll(Neuron,tr) = sum(ll);
        
    end
end
end

