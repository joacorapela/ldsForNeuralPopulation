
% This script opens summarymatfiles, breaks into train and test datasets 
% (to prepare for Leave-one-neuron out in trial-based fits)
% and saves as LONO-timestamp-matfilename in the same preprocessing folder 
% to double check: Laser indexing in assignToFold

clear all
%close all

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

%rootdir  = '/mnt/data/Mitra/figs/';
rootdir = '/mnt/javadzam/winstor/swc/mrsic_flogel/public/projects/MiJa_20160601_VisualLongRangeConnectivity/Ephys/figs/';

nFolds = 5;
Seed = 1;
type = 'exGrooming';
% options: 'all', 'onlyCorrect', 'exGrooming','onlyCorrect_exGrooming'
noTestSet = 0; % default 0

for animali = 8:length(animallist)
    animalname = animallist{animali};
    
    summarymatfile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/task_*.mat',animallist{animali},preprocessinglist{animali})));
    res = load(fullfile(summarymatfile.folder,summarymatfile.name));
    % !! IMPORTANT: make sure the correct trials are correct, given
    % that they can be from 2 different files.
    perffile = ...
    dir(fullfile(rootdir,sprintf('%s/preprocessing/%s/rev_perf*.mat',animallist{animali},preprocessinglist{animali})));
    if ~isempty(perffile)
        p = load(fullfile(summarymatfile.folder,perffile.name));
        res.correctgotrialind = p.correctgotrialind;
        res.correctnogotrialind = p.correctnogotrialind;
    end
    
    res.correcttrialind = [res.correctgotrialind,res.correctnogotrialind];
 
    
    % choose valid trials
    nAllTr = size(res.PAllOn,1);
    if strcmp(type,'onlyCorrect_exGrooming')
        validTrIdx = intersect(res.correcttrialind,res.nogroomingind);
    elseif strcmp(type,'onlyCorrect')
        validTrIdx = res.correcttrialind;
    elseif strcmp(type,'exGrooming')
        validTrIdx = res.nogroomingind;
    elseif strcmp(type,'all')
        validTrIdx = 1:nAllTr;
    end
    nTr = numel(validTrIdx);
    % define test and train trial sets
    % nTr = size(res.PAllOn,1); % the number of all trials  
    % if exclude grooming, better be done before next line (or can be done
    % later, but number of test/train trials within each fold will be more
    % variable  
    if ~noTestSet
        Fold = cell(1,nFolds);
        rng(Seed,'twister')
        indices = crossvalind('Kfold',nTr,nFolds);
        for FoldN = 1:nFolds
            Fold{FoldN}.testInd = validTrIdx(find(indices == FoldN))';
            Fold{FoldN}.trainInd = validTrIdx(find(~(indices == FoldN)))';
        end
    else
        for FoldN = 1:nFolds            
            Fold{FoldN}.trainInd = validTrIdx';
        end
    end
    resultsFilename = fullfile(summarymatfile.folder,...
        ['LONO_',type,'_',datestr(now,'yy_mm_dd_HH_MM_SS'),'_',summarymatfile.name]);    
    save(resultsFilename,'Fold','nFolds','Seed','type','noTestSet');

end
