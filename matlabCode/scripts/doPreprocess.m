
% saves v6 mat file and ini files for multiple animals in data/animalname

binsizesecs = 0.2;
laserduration = 0.15;

%
addpath('/mnt/data/Mitra/cache/repos/ini2struct')
addpath('/mnt/data/Mitra/cache/repos/struct2ini')
addpath('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/Rcode/scripts')


cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')

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

rootdir  = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/data/';
resultdir = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/';

for i = 1%:length(animallist)
    if ~isdir(fullfile(rootdir,animallist{i}))
        mkdir(fullfile(rootdir,animallist{i}))
    end
     if ~isdir(fullfile(resultdir,animallist{i},[num2str(1000*binsizesecs),'ms']))
        mkdir(fullfile(resultdir,animallist{i},[num2str(1000*binsizesecs),'ms']))
     end

    
    matfilename = ...
    dir(sprintf('/mnt/data/Mitra/figs/%s/preprocessing/%s/task_*.mat',animallist{i},preprocessinglist{i}));
    
    %%% v6 mat file
    saveDataSubset(animallist{i}, matfilename.name(1:end-4), preprocessinglist{i})
    
    FForFB = exptype{i};
    save(sprintf('../../data/%s/exptype.mat', animallist{i}),'FForFB');
    
    %%% ini file
    copyfile('../../data/exampleMouse/binLDStimeSeries.ini',...
        sprintf('../../data/%s/binLDStimeSeries.ini',animallist{i}))
    
    config = ini2struct(sprintf('../../data/%s/binLDStimeSeries.ini',animallist{i}));    
    config.filenames.matlabdatafilename = strrep(config.filenames.matlabdatafilename,'exampleMouse',animallist{i});
    config.filenames.matlabdatafilename = strrep(config.filenames.matlabdatafilename,'###',matfilename.name(1:end-4));
    config.filenames.savefilename = strrep(config.filenames.savefilename,'exampleMouse',...
        [animallist{i},'/',num2str(1000*binsizesecs),'ms']);
    config.filenames.savefilename = strrep(config.filenames.savefilename,'###',matfilename.name(1:end-4));
    
    % set binsize
    config.config_params.binsizesecs = binsizesecs;
    config.config_params.laserduration = laserduration;
 
    
    struct2ini(sprintf('../../data/%s/binLDStimeSeries.ini',animallist{i}),config)
    
    %%% bin
    animalname = animallist{i};
    run('doSaveBinnedLDSTimeSeries.m')
end
