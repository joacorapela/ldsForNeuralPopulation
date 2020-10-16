load task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_full.mat

V1Shaft1 = V1{1};
V1Shaft2 = V1{2};
LMShaft1 = LM{1};
LMShaft2 = LM{2};

save('task_2019-02-06_21-36-35_preprocessing_2019_04_05_15_04_39_ks2_subset_V6.mat', '-V6', 'V1Shaft1', 'V1Shaft2', 'LMShaft1', 'LMShaft2', 'nogotrialind', 'gotrialind', 'LaserDelay', 'PStepTimeStampOn', 'PStepTimeStampOff')
