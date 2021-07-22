clear all
%close all
cd('/mnt/data/Mitra/cache/repos/ldsForNeuralPopulation/matlabCode/scripts')


animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};

matfilenamelist = {'task_remade_from_stimlaser_perf_timeSeries.mat','task_remade_from_stimlaser_perf_timeSeries.mat',...
    'task_remade_from_stimlaser_perf_timeSeries.mat','task_remade_from_stimlaser_perf_timeSeries.mat',...
    'task_2019-12-18_21-53-06_preprocessing_2020_03_02_19_35_12_ks2_full_timeSeries.mat','task_2019-12-23_21-59-43_preprocessing_2020_03_02_19_58_50_ks2_full_timeSeries.mat',...
    'task_2019-12-27_21-36-59_preprocessing_2020_03_02_20_32_35_ks2_full_timeSeries.mat','task_2019-01-30_21-11-58_preprocessing_2020_03_02_13_10_37_ks2_full_timeSeries.mat',...
    'task_2019-02-06_21-36-35_preprocessing_2020_03_02_13_27_45_ks2_full_timeSeries.mat','task_remade_from_stimlaser_perf_timeSeries.mat',...
    'task_remade_from_stimlaser_perf_timeSeries.mat','task_remade_from_stimlaser_perf_timeSeries.mat',...
    'task_remade_from_stimlaser_perf_timeSeries.mat','task_2019-12-15_23-15-22_preprocessing_2020_03_02_14_31_53_ks2_full_timeSeries.mat'};



for animali = 1:length(animallist) % animal5 error
    animalname = animallist{animali};
    matfilename = matfilenamelist{animali};
    binsize = '20ms/';
    skipmsg = 1;
    nrep = 5;
    
    area = 'V1';
    
    try
        run('doAnalyze.m')
    catch
    end
    
    area = 'LM';
    try
        run('doAnalyze.m')
    catch
    end
    
end
