% addpath('/mnt/data/Mitra/cache/repos/ini2struct')

processAll(animalname)

function spikeCounts = getSpikeCounts(spikeSamples, breaks)
    nUnits = length(spikeSamples);
    spikeCounts = NaN(nUnits, length(breaks)-1);
    for n=1:nUnits
        spikeCounts(n,:) = histcounts(spikeSamples{n}, breaks);
    end
end

function binnedStimulus = getBinnedStimulus(stimOnSamples, stimOffSamples, breaks)
    binnedStimulus = zeros(length(breaks)-1, 1);
    for i=1:length(stimOnSamples)
        stimulatedBins = find(stimOnSamples(i)<=breaks & breaks<stimOffSamples(i));
        binnedStimulus(stimulatedBins) = 1.0;
    end
end

function processAll(animalname)
    iniFilename = sprintf('../../data/%s/binLDStimeSeries.ini',animalname);

    config = ini2struct(iniFilename);

    sRate = str2num(config.config_params.srate);
    binSizeSecs = str2num(config.config_params.binsizesecs);
    laserDuration = str2num(config.config_params.laserduration);
    matlabDataFilename = config.filenames.matlabdatafilename;
    saveFilename = config.filenames.savefilename;

    loadRes = load(matlabDataFilename);

    maxSpikeSample = 0;

    aux  = loadRes.V1Shaft1.SingleUnitSpikeTimes;
    nUnits = length(aux);
    v1Shaft1SpikeSamples = {};
    for n=1:nUnits
        v1Shaft1SpikeSamples{n} = aux{n};
        maxSpikeSample = max(maxSpikeSample, max(v1Shaft1SpikeSamples{n}));
    end

    aux  = loadRes.V1Shaft2.SingleUnitSpikeTimes;
    nUnits = length(aux);
    v1Shaft2SpikeSamples = {};
    for n=1:nUnits
        v1Shaft2SpikeSamples{n} = aux{n};
        maxSpikeSample = max(maxSpikeSample, max(v1Shaft2SpikeSamples{n}));
    end

    aux  = loadRes.LMShaft1.SingleUnitSpikeTimes;
    nUnits = length(aux);
    lmShaft1SpikeSamples = {};
    for n=1:nUnits
        lmShaft1SpikeSamples{n} = aux{n};
        maxSpikeSample = max(maxSpikeSample, max(lmShaft1SpikeSamples{n}));
    end

    aux  = loadRes.LMShaft2.SingleUnitSpikeTimes;
    nUnits = length(aux);
    lmShaft2SpikeSamples = {};
    for n=1:nUnits
        lmShaft2SpikeSamples{n} = aux{n};
        maxSpikeSample = max(maxSpikeSample, max(lmShaft2SpikeSamples{n}));
    end

    binSizeSamples = round(binSizeSecs*sRate);
    breaks = 0:binSizeSamples:(maxSpikeSample+binSizeSamples);
    v1Shaft1SpikeCounts = getSpikeCounts(v1Shaft1SpikeSamples, breaks);
    v1Shaft2SpikeCounts = getSpikeCounts(v1Shaft2SpikeSamples, breaks);
    lmShaft1SpikeCounts = getSpikeCounts(lmShaft1SpikeSamples, breaks);
    lmShaft2SpikeCounts = getSpikeCounts(lmShaft2SpikeSamples, breaks);

    stimOnSamples = loadRes.PStepTimeStampOn;
    stimOffSamples = loadRes.PStepTimeStampOff;
    goTrialIndices = loadRes.gotrialind;
    nogoTrialIndices = loadRes.nogotrialind;
    laserDelays = loadRes.LaserDelay/1000; %laser delays are in ms and I want them in sec

    goStimBinned = getBinnedStimulus(stimOnSamples(goTrialIndices), stimOffSamples(goTrialIndices), breaks);
    nogoStimBinned = getBinnedStimulus(stimOnSamples(nogoTrialIndices), stimOffSamples(nogoTrialIndices), breaks);
    laserOnsetSamples = stimOnSamples+round(laserDelays*sRate);
    laserOnsetSamples = laserOnsetSamples(find(~isnan(laserOnsetSamples)));
    laserStimBinned = getBinnedStimulus(laserOnsetSamples, laserOnsetSamples+round(laserDuration*sRate), breaks);

    timeSeries = struct(...
        'sRate', 1.0/binSizeSecs,...
        'breaks', breaks,...
        'v1Shaft1SpikeCounts', v1Shaft1SpikeCounts,...
        'v1Shaft2SpikeCounts', v1Shaft2SpikeCounts,...
        'lmShaft1SpikeCounts', lmShaft1SpikeCounts,...
        'lmShaft2SpikeCounts', lmShaft2SpikeCounts,...
        'goStim', goStimBinned,...
        'nogoStim', nogoStimBinned,...
        'laserStim', laserStimBinned,...
        'goStimOnSecs', stimOnSamples(goTrialIndices)/sRate,...
        'goStimOffSecs', stimOffSamples(goTrialIndices)/sRate,...
        'nogoStimOnSecs', stimOnSamples(nogoTrialIndices)/sRate,...
        'nogoStimOffSecs', stimOffSamples(nogoTrialIndices)/sRate,...
        'laserStimOnSecs', laserOnsetSamples/sRate,...
        'laserStimOffSecs', laserOnsetSamples/sRate+laserDuration);
    save(saveFilename, 'timeSeries');

end

