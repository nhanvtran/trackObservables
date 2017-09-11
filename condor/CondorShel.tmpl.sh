#!/bin/tcsh -f

setenv DISPLAY 0
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog:FASTJETLOC/bin

cd CMSSWBASE/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -csh`
rehash

xrdcp ANASUBLOC ${_CONDOR_SCRATCH_DIR}/anaSubstructure
xrdcp INDIR/FILE.dat ${_CONDOR_SCRATCH_DIR}/
cd ${_CONDOR_SCRATCH_DIR}
echo `grep 'HepMC' FILE.dat | wc` 

echo 'here'

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/uscms_data/d3/cvernier/TrackObservables/pheno/trackObservables/processing/2.06.09/lib

chmod 777 anaSubstructure
./anaSubstructure FILE.dat ./ 0 100000 CFG 
ls
xrdcp processed-FILE.dat-CFG.root root://cmseos.fnal.gov:///store/user/cvernier/OUTDIRFOLD/
rm *.dat *.root anaSubstructure 
