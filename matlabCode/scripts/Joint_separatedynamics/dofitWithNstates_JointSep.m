function [params ,seq ,varBound ,EStepTimes ,MStepTimes] = dofitWithNstates_JointSep(xDim,seq,Inference_handle,config,N_V1,N_LM)

uDim    = size(seq(1).u,1);
% xDim    = 9;
yDim    = size(seq(1).y,1);
T       = size(seq(1).T);
Trials  = length(seq);
maxIter = 100;
doff    = 0.0;

% if config.baselineU
%     for i =1:length(seq)
%         seq(i).u(3:end,:) = ~seq(i).u(3:end,:);
%     end
% end

fprintf('Max spike count:    %i \n', max(vec([seq.y])))
fprintf('Mean spike count:   %d \n', mean(vec([seq.y])))
fprintf('Freq non-zero bin:  %d \n', mean(vec([seq.y])>0.5))

%%% fit model

rng(config.RandSeed,'twister')

params = [];
% important: set flag to use external input
if uDim>0;params.model.notes.useB = true;end

%params = PLDSInitialize(seq, xDim, 'NucNormMin', params);
params = PLDSInitialize(seq, xDim, 'params', params);
% fprintf('Initial subspace angle:  %d \n', subspace(tp.model.C,params.model.C))

% params.seed = RandSeed;

%%% try giving C a mask
params.model.notes.useCMask = 1;
params.model.CMask = ones(size(params.model.C));
for rowNeuron = 1:N_V1
    params.model.CMask(rowNeuron,size(params.model.C,2)/2+1:end) = 0;
end
for rowNeuron = N_V1+1:N_V1+N_LM
    params.model.CMask(rowNeuron,1:size(params.model.C,2)/2) = 0;
end
% for rowNeuron = 1:N_V1
%     params.model.CMask(rowNeuron,4:6) = 0;
% end
% for rowNeuron = N_V1+1:N_V1+N_LM
%     params.model.CMask(rowNeuron,1:3) = 0;
% end

%params.opts.algorithmic.EMIterations = touchField(params.opts.algorithmic.EMIterations,'abortDecresingVarBound',false);

params.model.C = params.model.C.*params.model.CMask;
params.model.notes.maskB = 0; 
% temp: for 10 states
% toAdd=  blkdiag(0.1+zeros(5,5),0.1+zeros(5,5))- 0.2*eye(10);
% params.model.A = params.model.A + toAdd;
% 

%params.model.inferenceHandle = @PLDSVariationalInference;
params.model.inferenceHandle = Inference_handle; %@PLDSLaplaceInference;
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;
tic; [params seq varBound EStepTimes MStepTimes] = PopSpikeEM(params,seq); toc
% fprintf('Final subspace angle:  %d \n', subspace(tp.model.C,params.model.C))