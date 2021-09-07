
% set this variables before running
%nStates = 9;
area = 'V1';
%animali = 1;

addpath('./Aux')
addpath('./Analysis')
%addpath('/mnt/data/Mitra/cache/repos/ini2struct')
%addpath('/mnt/data/Mitra/cache/repos/struct2ini')
%addpath('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/Rcode/scripts')
%addpath('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

%rootdir  = '/mnt/data/Mitra/figs/';
%resultdir = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/';
rootdir  = '/nfs/winstor/mrsic_flogel/public/projects/MiJa_20160601_VisualLongRangeConnectivity/Ephys/figs/';
resultdir = '../../results/';
skipifexists = 1;
splitDelays = 0;
%cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};

preprocessinglist = {'2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48',...
    '2020_03_02_19_35_12','2020_03_02_19_58_50','2020_03_02_20_32_35',...
    '2020_03_02_13_10_37','2020_03_02_13_27_45',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19','2020_03_02_14_31_53'};

exptype = {'FF','FF','FF','FF',...
    'FF','FF','FF',...
    'FB','FB',...
    'FB','FB','FB','FB','FB'};

animalname = animallist{animali};
binSizems = 50;
binWinms = nan;% example: [500,1000]; pre-post stimulus window - nan: default [-1000 to 1000], custom windows not implemented yet
skipmsg = 1;

LONO.do = 1; % if 1 uses the train set only
LONO.fold = 1; % fold number to use for now keep a number, (not implemented) if 'all' uses and saves all: put an if below?

summarymatfile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/task_*.mat',animallist{animali},preprocessinglist{animali})));

if ~LONO.do
    savedir = fullfile(resultdir,animalname,'trial_based',[num2str(binSizems),'ms','Bins']);
    savename = [area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
else
    LONO.file = ...
        dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/LONO_*.mat',animallist{animali},preprocessinglist{animali})));
    % default: latest one
    if length(LONO.file) > 1
        LONO.file = LONO.file(end);
    end
    % add a summarymatfile name option for LONO - also save in lono folder,
    % name fold1
    savedir = fullfile(resultdir,animalname,'trial_based_LONO_MS',[num2str(binSizems),'ms','Bins']);
    savename = ['Fold',num2str(LONO.fold),'_',area,'_PLDSfitRes_',datestr(now,'yy_mm_dd_HH_MM_SS')];
    savename = ['nStates',num2str(nStates),'-',savename];   
end

if ~isdir(savedir)
    mkdir(savedir)
end
resultsFilename = [savedir,'/',savename,'.mat'];
resultsFigname = [savedir,'/',savename,'.pdf'];

if skipifexists
    if numel(dir([savedir,'/','*nStates',num2str(nStates),'*.mat']))
        return
    end
end


seq = buildTrialBasedSeq(summarymatfile, binSizems,binWinms,area,LONO,splitDelays);
config.binSizems= binSizems;
config.binWinms = binWinms;
config.area = area;
config.animalname = animalname;
config.summarymatfile = summarymatfile;
config.splitDelays = splitDelays;


codeRoot = '../../../pop_spike_dyn';
oldFolder = cd(codeRoot);
set_path
cd(oldFolder)

dbstop if error


clear resTr
FittedFold = LONO.fold;
try % it might error with some nsts
    [resTr.params ,resTr.seq ,resTr.varBound ,resTr.EStepTimes ,resTr.MStepTimes] = dofitWithNstates(nStates,seq);
    resTr.LONO = LONO;
    resTr.config = config;
    doLONO_and_ValLogLik;
    resTr.lono_trial_ll = trial_ll;
    resTr.lono_Avtrial_ll = mean(mean(trial_ll,2));
    resTr.test_trial_ll = test_trial_ll;
    resTr.test_Avtrial_ll = mean(mean(test_trial_ll,2));
catch
end
resTr.nStates = nStates;
% save trial_ll along with ns and some state variables
% also add option for repeating folds
ModelSelection = resTr;

save(resultsFilename,'ModelSelection')
cd(oldFolder)

