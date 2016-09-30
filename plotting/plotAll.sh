#!/bin/sh

cd /uscms_data/d3/ecoleman/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd /uscms_data/d3/ecoleman/TrackObservablesStudy/trackObservables/plotting/

anasubs=(nores_nogran res_nogran res_gran0p05 res_gran0p005 resx5_nogran resx100_nogran)

filehandles=(
processed-pythia82-fcc100-WW-pt5-50k 
processed-pythia82-lhc13-WW-pt1-50k 
processed-pythia82-fcc100-ZZ-pt5-50k 
processed-pythia82-lhc13-ZZ-pt1-50k 
processed-pythia82-fcc100-gg-pt5-50k 
processed-pythia82-lhc13-gg-pt1-50k 
processed-pythia82-fcc100-tt-pt5-50k 
processed-pythia82-lhc13-tt-pt1-50k 
processed-pythia82-fcc100-qq-pt5-50k 
processed-pythia82-lhc13-qq-pt1-50k)

mkdir plots
for ana in ${anasubs[*]} ; do
    mkdir plots/${ana}
    rm plots/${ana}/*
for fileh in ${filehandles[*]} ; do
    echo " "
    echo " "
    echo " "
    echo "Working on ${fileh}-${ana}.root"
    echo " "
    echo
    nohup python quickPlotter.py -b --basedir ../processing --ana ${ana} -o ./plots/${ana}/ > ${ana}.txt &
done
done

for ana in ${anasubs[*]} ; do
    scp plots/${ana}/* ecoleman@lxplus.cern.ch:~/www/TrackObservablesStudy/SmearPlots/09_02_2016/${ana}/
done
