#!/bin/tcsh -f

# modify path and setup basic env
setenv DISPLAY 0
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog:FASTJETLOC/bin
source /cvmfs/cms.cern.ch/cmsset_default.csh

# setup cmssw environment
xrdcp CMSSWLOC/CMSSWVER.tar.gz ${_CONDOR_SCRATCH_DIR}/cmssw.tar.gz
tar -xf cmssw.tar.gz
rm -rf cmssw.tar.gz
cd ${_CONDOR_SCRATCH_DIR}/CMSSWVER/src/
eval `scramv1 runtime -csh`

# copy analysis files
xrdcp ANASUBLOC/anaSubstructure ${_CONDOR_SCRATCH_DIR}/anaSubstructure
xrdcp INDIR/FILE.lhe ${_CONDOR_SCRATCH_DIR}/
xrdcp ANASUBLOC/lhc14-pythia8-4C-minbias-nev100.pu14.gz ${_CONDOR_SCRATCH_DIR}/

# run analysis
cd ${_CONDOR_SCRATCH_DIR}
chmod 777 anaSubstructure
./anaSubstructure FILE ./ MINEV MAXEV CFG TAG
xrdcp processed-FILE-TAG.root root://cmseos.fnal.gov:///store/user/USER/OUTDIRFOLD/
rm -rf CMSSWVER/
rm *.lhe *.root *.pu14.gz anaSubstructure 
