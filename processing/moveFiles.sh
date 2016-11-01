#!/bin/sh

cd /uscms_data/d3/ecoleman/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd /uscms_data/d3/ecoleman/TrackObservablesStudy/trackObservables/processing/

for fileh in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/ntran/trackObservables/samples/ | grep lhe) ; do
    echo " "
    echo "Working on ${fileh}"
    xrdcp ${fileh} root://cmseos.fnal.gov:///store/user/ecoleman/TrackObservablesStudy/fromMarat/
done

