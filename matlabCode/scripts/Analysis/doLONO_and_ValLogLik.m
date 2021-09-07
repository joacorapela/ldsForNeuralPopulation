% uncomment for a single example:
%cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/VL61/trial_based_LONO/50msBins')
%Training = 'Fold1_V1_PLDSfitRes_21_08_10_14_28_34.mat';
%resTr=load(Training);
%FittedFold = str2num(strtok(strtok(Training,'Fold'),'_'));
doPlot = 0;

% corresponding LONO params file:
LONOparams = load(fullfile(resTr.LONO.file.folder,resTr.LONO.file.name));
% corresponding LONO params file:

% this is the fold used for training
%%% IMPORTANT: output this to double check using the correct train/test set
% or try shifting and see the improvement in accuract of test


nNeurons = size(resTr.params.model.C,1);
numTestTrials = length(LONOparams.Fold{FittedFold}.testInd);

trial_ll = nan(nNeurons,numTestTrials);
% if doSplit=1, there could be less trials than expected (the ones with delays 
% outside of discrete range are ignored. 
for heldoutN = 1:nNeurons % id of held out neuron
    tll = heldout_loglike(resTr,LONOparams,FittedFold,heldoutN,doPlot);
    trial_ll(heldoutN,1:size(tll,2)) = tll;
end
% removing nan trials: from the loop above if there are less trials than
% expected 
trial_ll = trial_ll(:,find(~isnan(sum(trial_ll,1))));
% ! careful that trial_ll values were summed over time bins, so, their
% value is sensitive to the number of time bins !
% averaged over trials and neuorns
mean(mean(trial_ll,2))

%%%%% loglikelihood of prediction on the test set (no heldout)


%test_trial_ll = nan(nNeurons,numTestTrials);
test_trial_ll = test_loglike(resTr,LONOparams,FittedFold,doPlot);
