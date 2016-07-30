#!/bin/tcsh -f

setenv DISPLAY 0
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog:FASTJETLOC/bin

cd CMSSWBASE/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -csh`
rehash

cp ANASUBLOC ${_CONDOR_SCRATCH_DIR}/
cd ${_CONDOR_SCRATCH_DIR}

./anaSubstructure FILE INDIR/ MINEV MAXEV TAG
