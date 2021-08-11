% cd to the right folder
Training = 'Fold1_V1_PLDSfitRes_21_08_10_14_28_34.mat';
resTr=load(Training);

% corresponding LONO params file:
LONOparams = load(fullfile(resTr.LONO.file.folder,resTr.LONO.file.name));
% corresponding LONO params file:

% this is the fold used for training
%%% IMPORTANT: output this to double check using the correct train/test set
% or try shifting and see the improvement in accuract of test
FittedFold = str2num(strtok(strtok(Training,'Fold'),'_'));

% make y of the test dataset, leaving out heldoutN
heldoutN = 4; % id of held out neuron
% seq_one will only get a y_orig field which is its real spiking rate
[seq_minusOne,seq_One] = buildTrialBasedSeq_heldout(resTr.config.summarymatfile, resTr.config.binSizems,resTr.config.binWinms,...
    resTr.config.area,LONOparams.Fold{FittedFold}.testInd,heldoutN);

% inference: get x - before here, prepare the loading matrix (C) and d : eliminate the corresponding row)
% check if anything else is needed
resTs_minusOne = resTr;
resTs_minusOne.params.model.C(heldoutN,:) = [];
resTs_minusOne.params.model.d(heldoutN) = [];
resTs_minusOne.seq = seq_minusOne;
[resTs_minusOne.seq,~] = resTr.params.model.inferenceHandle(resTs_minusOne.params,resTs_minusOne.seq); % check what it wants from seq 
%                                                                                  % (it has now a field named yorig with one row, no u)                                                                                  % (it has now a field named yorig with one row, no u) 


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


figure;plot(nanmean(pred,1));hold on;plot(nanmean(yOrig,1));legend({'prediction','original trace'})
% why the offset?

% measure the difference between pred and yOrig

