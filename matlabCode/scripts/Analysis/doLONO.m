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
for heldoutN = 1:nNeurons % id of held out neuron
    trial_ll(heldoutN,:) = heldout_loglike(resTr,LONOparams,FittedFold,heldoutN,doPlot);
end
% ! careful that trial_ll values were summed over time bins, so, their
% value is sensitive to the number of time bins !
% averaged over trials and neuorns
mean(mean(trial_ll,2))

