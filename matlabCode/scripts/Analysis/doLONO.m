% cd to the right folder
Training = 'Fold1_V1_PLDSfitRes_21_08_09_14_27_24.mat';
resTr=load(Training);

% corresponding LONO params file:
LONOparams = load(fullfile(resTr.LONO.file.folder,resTr.LONO.file.name));
% corresponding LONO params file:

% this is the fold used for training
%%% IMPORTANT: output this to double check using the correct train/test set
% or try shifting and see the improvement in accuract of test
FittedFold = str2num(strtok(strtok(Training,'Fold'),'_'));

% make y of the test dataset, leaving out heldoutN
heldoutN = 3; % id of held out neuron
% seq_one will only get a y_orig field which is its real spiking rate
[seq_minusOne,seq_One] = buildTrialBasedSeq_heldout(resTr.config.summarymatfile, resTr.config.binSizems,resTr.config.binWinms,...
    resTr.config.area,LONOparams.Fold{FittedFold}.testInd,heldoutN);

% inference: get x - before here, prepare the loading matrix (C) and d : eliminate the corresponding row)
% check if anything else is needed
resTr_minusOne = resTr;
resTr_minusOne.params.model.C(heldoutN,:) = [];
resTr_minusOne.params.model.d(heldoutN) = [];
[seq_One,~] = resTr.params.model.inferenceHandle(resTr_minusOne.params,seq_minusOne);

% for each trial: compare seq_One(1).y and the actual neural trace
% seq_One(1).yOrig

