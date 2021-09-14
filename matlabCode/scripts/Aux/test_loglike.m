
function trial_ll = test_loglike(resTr,LONOparams,FittedFold,doPlot)

% make y of the test dataset
seq = buildTrialBasedSeq_test(resTr.config.summarymatfile, resTr.config.binSizems,resTr.config.binWinms,...
    resTr.config.area,LONOparams.Fold{FittedFold}.testInd,resTr.config.splitDelays);

% inference: get x
resTs = resTr;
resTs.seq = seq;

filtered = 0; % if 1 filtered states, if not smoothed states. 

if filtered
    x_filt = zeros(1,40,resTr.nStates);
    for i =10:length(x_filt) % start from 10 somehow
        temp = resTs.seq;
        for tr = 1:length(resTs.seq)
            temp(tr).y = resTs.seq(tr).y(:,1:i);
            temp(tr).u = resTs.seq(tr).u(:,1:i);
            temp(tr).T= i;
        end
        [temp,~] = resTr.params.model.inferenceHandle(resTs.params,temp);
        
        for tr = 1:length(resTs.seq)
            x_filt(tr,i,:) = temp(tr).posterior.xsm(:,i);
        end
    end
    for tr = 1:length(resTs.seq)
        resTs.seq(tr).posterior.xsm = squeeze(x_filt(tr,:,:))';
    end
    
else
    [resTs.seq,~] = resTr.params.model.inferenceHandle(resTs.params,resTs.seq);
end


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
        if filtered
            lambda = pred(tr,10:40);
            x = yOrig(tr,10:40);
            n = length(yOrig(tr,10:40));
        else
            lambda = pred(tr,:);
            x = yOrig(tr,:);
            n = length(yOrig(tr,:));
        end
        
        % this is a 1*n array, each element loglikelihood of one timepoint
        ll = -lambda + x.*log(lambda) - log(factorial(x));
        % ll of individual time bins are summed to make the total ll of the
        % trial
        trial_ll(Neuron,tr) = sum(ll);
        
    end
end
end

