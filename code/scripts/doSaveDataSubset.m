load ../../data/VL61/stimlaserL_2018-01-12_18-17-46_preprocessing_2018_03_16_16_53_13.mat

if(iscell(nogotrialind))
    nogotrialind = nogotrialind{1};
end

if(iscell(gotrialind))
    gotrialind = gotrialind{1};
end

if(iscell(LaserDelay))
    LaserDelay = LaserDelay{1};
end

if(iscell(PStepTimeStampOn))
    PStepTimeStampOn = PStepTimeStampOn{1};
end

if(iscell(PStepTimeStampOff))
    PStepTimeStampOff = PStepTimeStampOff{1};
end

V1Shaft1 = V1{1};
V1Shaft2 = V1{2};
LMShaft1 = LM{1};
LMShaft2 = LM{2};

save('stimlaserL_2018-01-12_18-17-46_preprocessing_2018_03_16_16_53_13_subset_V6.mat', '-v6', 'V1Shaft1', 'V1Shaft2', 'LMShaft1', 'LMShaft2', 'nogotrialind', 'gotrialind', 'LaserDelay', 'PStepTimeStampOn', 'PStepTimeStampOff')
% save stimlaserL_2018-01-12_18-17-46_preprocessing_2018_03_16_16_53_13_subset_V6.mat -v6 V1Shaft1 V1Shaft2 LMShaft1 LMShaft2 nogotrialind gotrialind LaserDelay PStepTimeStampOn PStepTimeStampOff
% save stimlaserL_2018-01-12_18-17-46_preprocessing_2018_03_16_16_53_13_subset.mat V1Shaft1 V1Shaft2 LMShaft1 LMShaft2 nogotrialind gotrialind LaserDelay PStepTimeStampOn PStepTimeStampOff
% save('stimlaserL_2018-01-12_18-17-46_preprocessing_2018_03_16_16_53_13_subset_V6.mat', '-v6', 'V1Shaft1')
