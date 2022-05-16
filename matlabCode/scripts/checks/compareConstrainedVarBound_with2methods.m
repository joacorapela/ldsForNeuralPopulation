%% run doFit_Allanimals_JointSepCon.m with a break point at the first line of  dofitWithNstates_JointSepCon.m, then run this code
%%
uDim    = size(seq(1).u,1);
yDim    = size(seq(1).y,1);
T       = size(seq(1).T);
Trials  = length(seq);
maxIter = 50;%20
doff    = 0.0;

%%
clear out_nc out_c_m1 out_c_m2
%% no B constraint
params = [];
params.model.notes.useB = true;
params.model.notes.maskB = false;
params.model.notes.maskBmethod = '';
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


params.model.inferenceHandle = @PLDSLaplaceInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_nc.params, out_nc.seq, out_nc.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%% B constraint m1
params = [];
params.model.notes.useB = true;
params.model.notes.maskB = true;
params.model.notes.maskBmethod = 'method1';

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


params.model.inferenceHandle = @PLDSLaplaceInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_c_m1.params, out_c_m1.seq, out_c_m1.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%% B constraint m2
params = [];
params.model.notes.useB = true;
params.model.notes.maskB = true;
params.model.notes.maskBmethod = 'method2';

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


params.model.inferenceHandle = @PLDSLaplaceInference; % changed from default
params.opts.algorithmic.EMIterations.maxIter     = maxIter;
params.opts.algorithmic.EMIterations.maxCPUTime  = inf;

tic; [out_c_m2.params, out_c_m2.seq, out_c_m2.varBound, ~, ~] = PopSpikeEM(params,seq); toc
%% add no b? meaninful?
%% 
% save
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/ConstraintTest');
save('res50_l_2m_q.mat','out_nc','out_c_m1','out_c_m2')
%%
figure;plot(out_nc.varBound,'k')
hold on;plot(out_c_m1.varBound,'r')
hold on;plot(out_c_m2.varBound,'c')

legend({'no constraint','constraint m1','constraint m2'})
%%
iter = 1;
climm = [-.1 .1];
figure;subplot(1,2,1);imagesc(out_c_m1.params.model.difA(:,:,iter));caxis(climm);subplot(1,2,2);imagesc(out_c_m1.params.model.difB(:,:,iter));caxis(climm);colormap('gray')
figure;subplot(1,2,1);imagesc(out_c_m2.params.model.difA(:,:,iter));caxis(climm);subplot(1,2,2);imagesc(out_c_m2.params.model.difB(:,:,iter));caxis(climm);colormap('gray')
%% 
figure;plot(squeeze(out_nc.params.model.Qdyn),'k')
hold on;plot(squeeze(out_c_m1.params.model.Qdyn),'r')
hold on;plot(squeeze(out_c_m2.params.model.Qdyn),'c')

legend({'no constraint','constraint m1','constraint m2'})