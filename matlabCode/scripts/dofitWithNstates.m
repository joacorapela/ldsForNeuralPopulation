function [params ,seq ,varBound ,EStepTimes ,MStepTimes] = dofitWithNstates(xDim,seq,Inference_handle,config)

uDim    = size(seq(1).u,1);
% xDim    = 9;
yDim    = size(seq(1).y,1);
T       = size(seq(1).T);
Trials  = length(seq);
maxIter = 100;
doff    = 0.0;

if config.baselineU
    if ~config.dynInput
        for i =1:length(seq)
            seq(i).u(3:end,:) = ~seq(i).u(3:end,:);
        end
    else
        for i =1:length(seq)
            seq(i).u((end-15):end,:) = ~seq(i).u((end-15):end,:);
        end
    end
end

fprintf('Max spike count:    %i \n', max(vec([seq.y])))
fprintf('Mean spike count:   %d \n', mean(vec([seq.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([seq.y])>0.5))
fprintf('init params:  %s %d\n', config.InitType,config.RandSeed)
%%% fit model

params = [];
% important: set flag to use external input
if uDim>0;params.model.notes.useB = true;end

%params = PLDSInitialize(seq, xDim, 'NucNormMin', params);


rng(config.RandSeed,'twister')
params = PLDSInitialize(seq, xDim, config.InitType, params);

% fprintf('Initial subspace angle:  %d \n', subspace(tp.model.C,params.model.C))
if config.AfromJoint
    params.model.notes.learnA = 0;
    Jdir = dir(['/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/',config.animalname,'/Joint_trial_based_split_0/50msBins/','*.mat']);
    J = load(fullfile(Jdir(end).folder,Jdir(end).name));
    % temp: last one is 16 states
    if size(J.params.model.A,1) == 16
        params.model.A = J.params.model.A(1:8,1:8);
    else
        error
    end
end

%params.model.inferenceHandle = @PLDSVariationalInference;
params.model.inferenceHandle = Inference_handle; %@PLDSLaplaceInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))