
% compares variance of the data along the dimensions of C (low dimensional
% latent space) to the variance of data along its principal components


cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results')

animalname = 'VL63';
binsize= '200ms'; % or '' for old files
timeSeriesFilename = 'task_remade_from_stimlaser_perf_timeSeries.mat';
area = 'V1';

trainDurSecs = 300;%180
analysisStartTimeSecs = (4)*trainDurSecs; % choose the model to load here (start time)

cd(fullfile(animalname,binsize))

% load timeseries file
timeSeries = load(timeSeriesFilename);
timeSeries = timeSeries.timeSeries;
sRate = timeSeries.sRate;
trainSamples = (analysisStartTimeSecs*sRate)+(1:(trainDurSecs*sRate));

if strcmp(area,'V1')
    spikeCounts = [timeSeries.v1Shaft1SpikeCounts;timeSeries.v1Shaft2SpikeCounts];
elseif strcmp(area,'LM')
    spikeCounts = [timeSeries.lmShaft1SpikeCounts;timeSeries.lmShaft2SpikeCounts];%timeSeries.v1Shaft1SpikeCounts;
end

spikeCounts = spikeCounts(:,trainSamples);

% load PLDS results
cd(num2str(trainDurSecs))
PLDSresFiles = dir([area,sprintf('*_%s.mat',num2str(analysisStartTimeSecs))]);
res = load(PLDSresFiles.name);


% size(res.params.model.C)
% size(pca(spikeCounts'))

[Vs,score,lambda] = pca(spikeCounts','Centered',1); % lambda = var(spikeCounts'*Vs)
% normalize columns of res.params.model.C (?)
for i = 1:size(res.params.model.C,2)
    res.params.model.C(:,i) = res.params.model.C(:,i)./norm(res.params.model.C(:,i));
end

[VarAlongC,~] = sort(var(spikeCounts'*res.params.model.C),'descend');
% take as many pcs as the number of states (size(res.params.model.C,2))
VarAlongPC = var(spikeCounts'*Vs(:,1:size(res.params.model.C,2)));

figure;plot(VarAlongC,'b')
hold on;plot(VarAlongPC,'k')

% %% project on pcs
% for i = 1:length(ord)
%     res.params.model.C(:,ord(i))'*Vs(:,i)
% end
