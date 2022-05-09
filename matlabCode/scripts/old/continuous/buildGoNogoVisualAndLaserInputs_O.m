
function allInputs = buildGoNogoVisualAndLaserInputs(goStim, nogoStim, laserStim, inputMemory)
    goStimInputs = buildInputsBlock(goStim, inputMemory);
    nogoStimInputs = buildInputsBlock(nogoStim, inputMemory);
    laserInputs = buildInputsBlock(laserStim, inputMemory);
    goLaserInputs = goStimInputs.*laserInputs;
    nogoLaserInputs = nogoStimInputs.*laserInputs;
    allInputs = [goStimInputs; nogoStimInputs; goLaserInputs; nogoLaserInputs];
end

function inputs = buildInputsBlock(stim, inputMemory)
    N = length(stim);
    inputs = nan(inputMemory+1,N);
    for i=0:(inputMemory)
        inputs(i+1,:) = [zeros(1,i), stim(1:(N-i))];
    end
end

