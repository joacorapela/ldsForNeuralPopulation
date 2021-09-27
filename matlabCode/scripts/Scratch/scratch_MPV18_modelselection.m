% resTr

cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/cluster/MPV18_2/trial_based_MS_split_1/17msBins')
%load('nStates9-Fold1_V1_PLDSfitRes_21_09_09_15_42_40.mat')
%resTr = ModelSelection;

doPlot = 0;
FittedFold = 1; 

%warning ('off','all');

% by default reads tghe latest version
ll = nan(1,50);
for st = 1:numel(test_ll)
    st
    File = dir(['nStates',num2str(st),'-','*.mat']);
    % if many, last one
    if numel(File) > 1
       File = File(end);
    end
    if numel(File)
        res = load(File.name);
        resTr = res.ModelSelection;
        resTr.config.summarymatfile.folder = '/mnt/data/Mitra/figs/MPV18_2/preprocessing/2020_03_02_13_27_45';
        resTr.LONO.file.folder =  '/mnt/javadzam/winstor/swc/mrsic_flogel/public/projects/MiJa_20160601_VisualLongRangeConnectivity/Ephys/figs/MPV18_2/preprocessing/2020_03_02_13_27_45';
        LONOparams = load(fullfile(resTr.LONO.file.folder,resTr.LONO.file.name));
        doLONO;
        ll(st) = nanmean(nanmean(trial_ll));
    end
end

figure;plot(1:25,ll,'.-');
