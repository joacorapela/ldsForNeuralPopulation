#!/bin/csh

set filename='../../data/exampleMouse/v1Shaft1_modelSelection_tmp.txt'
set analysisStartTimeSecs=(0  180  360  540  720  900 1080 1260 1440 1620 1800 1980 2160 2340 2520 2700 2880 3060 3240 3420 3600 3780 3960)
set trainDurSecs=(180)
set validationDurSecs=(60)
set stateDims=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)
set initMethods=(PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA PPCA)
# set obsMemories=(0.4)
set obsMemories=(0.0)

foreach analysisStartTime ($analysisStartTimeSecs)
    foreach trainDur ($trainDurSecs)
        foreach validationDuration ($validationDurSecs)
            set i = 1
            while ($i <= ${#stateDims})
                foreach stateMemory ($stateMemories)
                    foreach obsMemory ($obsMemories)
                        echo $analysisStartTime $trainDur $validationDuration $stateDims[$i] $stateMemory $obsMemory $initMethods[$i] >> $filename
                    end
                end
                @ i++
            end
        end
    end
end
