#!/bin/sh

cd /uscms_data/d3/${USER}/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd /uscms_data/d3/${USER}/TrackObservablesStudy/trackObservables/processing/

anasubs=(nores_nogran res_nogran resx5_nogran resx100_nogran res_gran0p05 res_gran0p005)
filehandles=(processed-pythia82-fcc100-WW-pt5-50k processed-pythia82-lhc13-WW-pt1-50k)

for ana in ${anasubs[*]} ; do
for fileh in ${filehandles[*]} ; do
    echo " "
    echo " "
    echo " "
    echo "Working on ${fileh}-${ana}.root"
    echo " "
    hadd ${fileh}-${ana}.root  `xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/${ana} | grep root | grep ${fileh}`
done
done


