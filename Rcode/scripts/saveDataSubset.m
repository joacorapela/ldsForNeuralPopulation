function saveDataSubset(mouseName, dataFilenameInfix, varargin)
    
    if nargin > 2
         inputFilename = ...
        ['/mnt/data/Mitra/figs/',mouseName,'/preprocessing/',varargin{1},'/',dataFilenameInfix];
    else
        inputFilename = sprintf('../../data/%s/%s.mat', mouseName, dataFilenameInfix);
    end
    
    outputFilename = sprintf('../../data/%s/%s_subset_V6.mat', mouseName, dataFilenameInfix);

    load(inputFilename)

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

    save(outputFilename, '-v6', 'V1Shaft1', 'V1Shaft2', 'LMShaft1', 'LMShaft2', 'nogotrialind', 'gotrialind', 'LaserDelay', 'PStepTimeStampOn', 'PStepTimeStampOff')
