
uDim    = size(seq(1).u,1);
yDim    = size(seq(1).y,1);
T       = size(seq(1).T);
Trials  = length(seq);
maxIter = 5;%20
doff    = 0.0;

%%
clear out_nc out_c out_nb out_c_uA out_c_uAQ
%% no B constraint
params = [];
params.model.notes.useB = true;
params = PLDSInitialize(seq, xDim, 'NucNormMin', params); % changed from default
%%% try giving C a mask
params.model.notes.useCMask = 1;
params.model.CMask = ones(size(params.model.C));
for rowNeuron = 1:N_V1
    params.model.CMask(rowNeuron,size(params.model.C,2)/2+1:end) = 0;
end
for rowNeuron = N_V1+1:N_V1+N_LM
    params.model.CMask(rowNeuron,1:size(params.model.C,2)/2) = 0;
end
params.model.C = params.model.C.*params.model.CMask;

params.model.notes.maskB = 0;

params.model.inferenceHandle = @PLDSVariationalInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_nc.params, out_nc.seq, out_nc.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%% B constraint default
params = [];
params.model.notes.useB = true;
params = PLDSInitialize(seq, xDim, 'NucNormMin', params); % changed from default
%%% try giving C a mask
params.model.notes.useCMask = 1;
params.model.CMask = ones(size(params.model.C));
for rowNeuron = 1:N_V1
    params.model.CMask(rowNeuron,size(params.model.C,2)/2+1:end) = 0;
end
for rowNeuron = N_V1+1:N_V1+N_LM
    params.model.CMask(rowNeuron,1:size(params.model.C,2)/2) = 0;
end
params.model.C = params.model.C.*params.model.CMask;

params.model.notes.maskB = 1;

params.model.inferenceHandle = @PLDSVariationalInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_c.params, out_c.seq, out_c.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%% B conatraint update_A
params = [];
params.model.notes.useB = true;
params = PLDSInitialize(seq, xDim, 'NucNormMin', params); % changed from default
%%% try giving C a mask
params.model.notes.useCMask = 1;
params.model.CMask = ones(size(params.model.C));
for rowNeuron = 1:N_V1
    params.model.CMask(rowNeuron,size(params.model.C,2)/2+1:end) = 0;
end
for rowNeuron = N_V1+1:N_V1+N_LM
    params.model.CMask(rowNeuron,1:size(params.model.C,2)/2) = 0;
end
params.model.C = params.model.C.*params.model.CMask;

params.model.notes.maskB = 1;
params.replaceA = 1;

params.model.inferenceHandle = @PLDSVariationalInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_c_uA.params, out_c_uA.seq, out_c_uA.varBound, ~, ~] = PopSpikeEM(params,seq); toc


%% Constraint update A and Q
params = [];
params.model.notes.useB = true;
params = PLDSInitialize(seq, xDim, 'NucNormMin', params); % changed from default
%%% try giving C a mask
params.model.notes.useCMask = 1;
params.model.CMask = ones(size(params.model.C));
for rowNeuron = 1:N_V1
    params.model.CMask(rowNeuron,size(params.model.C,2)/2+1:end) = 0;
end
for rowNeuron = N_V1+1:N_V1+N_LM
    params.model.CMask(rowNeuron,1:size(params.model.C,2)/2) = 0;
end
params.model.C = params.model.C.*params.model.CMask;

params.model.notes.maskB = 1;
params.replaceA = 2;

params.model.inferenceHandle = @PLDSVariationalInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_c_uAQ.params, out_c_uAQ.seq, out_c_uAQ.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%%
% save
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/ConstraintTest');
save('res5_v.mat','out_nc','out_c','out_c_uAQ','out_c_uA')
%%
figure;plot(out_nc.varBound,'k')
hold on;plot(out_c.varBound,'r')
hold on;plot(out_c_uA.varBound,'c')
hold on;plot(out_c_uAQ.varBound,'g')

legend({'no constraint','constraint','uA','uAQ'})
