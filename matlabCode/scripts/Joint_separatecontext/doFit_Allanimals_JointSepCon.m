% Fits and saves the models for all animals.

% clear all
%close all


cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')
addpath('./Aux')
addpath('./Analysis')



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

rootdir  = '/mnt/data/Mitra/figs/';
resultdir = '/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/results/';

nStates = 16; % 1:35 for model selection
%splitDelays = 0; % default 0
Inference_handle = @PLDSLaplaceInference;
binSizems = 10; % default 50ms
binWinms = nan;% example: [500,1000]; pre-post stimulus window - nan: default [-1000 to 1000], custom windows not implemented yet
skipmsg = 1;
doSavefig = 1;
doSaveres = 1;
baselineU = 0; % 0 - option for 1 not implemented
% trialType = 'onlyCorrect_exGrooming_go';%'onlyCorrect_exGrooming'; % leave empty if all trials.


LONO.do = 0; % 1 option not implemented for joint yet
LONO.fold = 1; % fold number to use for now keep a number, (not implemented) if 'all' uses and saves all: put an if below?

% default only fb: changed to 1 to go through all
% However, for reversed stimuli (long range influnce), should separate
% animals (only indirect area)
for RandSeed =0:10
    for animali = 8%:13%length(animallist)
        animalname = animallist{animali};
        close all
        %try
        trialType = 'onlyCorrect_exGrooming_go';
        run('doFit_JointSepCon.m')
        clearvars -except animallist preprocessinglist exptype rootdir resultdir nStates Inference_handle binSizems binWinms skipmsg doSavefig doSaveres baselineU LONO RandSeed animali animalname trialType
        
        
        trialType = 'onlyCorrect_exGrooming_nogo';
        run('doFit_JointSepCon.m')
        clearvars -except animallist preprocessinglist exptype rootdir resultdir nStates Inference_handle binSizems binWinms skipmsg doSavefig doSaveres baselineU LONO RandSeed animali animalname trialType

        %catch
        %end       
        
    end
end

