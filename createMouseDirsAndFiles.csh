#!/bin/csh

if ($#argv != 1) then
    echo "Usage $0 mouseName"
    goto done
endif

set mouseName = $1
 
cp -r figures/exampleMouse figures/{$mouseName}
cp -r log/exampleMouse log/{$mouseName}
cp -r results/exampleMouse results/{$mouseName}
cp -r slurmOutputs/exampleMouse slurmOutputs/{$mouseName}

cp -r code/scripts/exampleMouse code/scripts/{$mouseName}
sed -i "s/exampleMouse/$mouseName/g" code/scripts/{$mouseName}/*.csh

cp -r data/exampleMouse data/{$mouseName}
sed -i "s/exampleMouse/$mouseName/g" data/{$mouseName}/*.ini data/{$mouseName}/*.txt

done:
 exit 0
