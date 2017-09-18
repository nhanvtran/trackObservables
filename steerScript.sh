#!/bin/bash

if [[ "$1" == "" ]]; then
    echo "*****************************************"
    echo "*                                       *"
    echo "* Track Observables Study: Steer Script *"
    echo "*                                       *"
    echo "*****************************************"
    echo "* MAKE   - build executable             *"
    echo "* JOBS   - launch condor jobs           *"
    echo "* MERGE  - merge outputs from EOS       *"
    echo "*****************************************"
    echo ""
    
    exit 1
fi

cd processing
procdir=$PWD

# Perfect
# H0.05
# H0.01
# F1 
# High-res
# CMS-like EB
# CMS-like EE
# Ultra-high res
# F2
#anasubs=(r0_h0_e0 
#r05_h05_e005 
#r05_h01_e005 
#r05_h01_e005_t500 
#r05_h002_e005_t500
# r1_h022_e050_t110
# r1_h022_e0175_t110
#r05_h005_e005_t500
#r05_h002_e001_t500) 

#anasubs=(r1_h022_e0175)
#anacfgs=(r1.0:h0.022:e0.0175)

#anasubs=(r1_h022_e0175_t220)
#anacfgs=(r1.0:h0.022:e0.0175:t220.0)

#anasubs=(r0_h0_e0)
#anacfgs=(z)

#anasubs=(r1_h022_e0175_t220_q10)
#anacfgs=(r1.0:h0.022:e0.0175:t220.0:q10)

# t220 F1 F2
anasubs=(r_h_e)
anacfgs=(r1.0:h0.022:t220.0)

#anasubs=(r1_h022_e0175_t500)
##r1_h022_e005_t110) 
#anacfgs=(r1:h0.022:e0.0175:t500.0)
##r1:h0.022:e0.005:t110.0) 

#anacfgs=(z 
#r0.5:h0.050:e0.005 
#r0.5:h0.010:e0.005
#r0.5:h0.010:e0.005:t500.0
#r0.5:h0.002:e0.005:t500.0
#r1.0:h0.022:e0.0500:t110.0  
#r1.0:h0.022:e0.0175:t110.0  
#r0.5:h0.005:e0.005:t500.0 
#r0.5:h0.002:e0.001:t500.0) 

filehandles=(
WJetsToQQ_HT-600toInf_tarball.tar.xz 
ZJetsToQQ_HT600toInf_gridpack.tar.gz
)

####################################################################
#             Nothing should change below this line                #
####################################################################


case $1 in

MAKE )
############################# MAKE #################################
echo "Making executables"

make
eval `eos root://cmseos.fnal.gov rm /store/user/cvernier/TrackObservablesStudy/anaSubstructure`
xrdcp anaSubstructure root://cmseos.fnal.gov:///store/user/cvernier/TrackObservablesStudy/
#xrdcp lhc14-pythia8-4C-minbias-nev100.pu14.gz root://cmseos.fnal.gov:///store/user/cvernier/TrackObservablesStudy/

#eval `eos root://eoscms.cern.ch rm /store/cmst3/group/exovv/precision/`
#${USER}/TrackObservablesStudy/anaSubstructure`
#xrdcp anaSubstructure root://eoscms.cern.ch:///store/group/phys_exotica/dijet/dazsle/TrackObservablesStudy/
#xrdcp lhc14-pythia8-4C-minbias-nev100.pu14.gz root://eoscms.cern.ch:///store/user/${USER}/TrackObservablesStudy/

;;

JOBS)
############################# JOBS #################################
echo "Running jobs"

i=0
for ana in ${anasubs[*]} ; do
    # create output directory, if it does not exist
    eval `eos root://cmseos.fnal.gov mkdir /store/user/cvernier/${ana}`
    #eval `eos root://eoscms.cern.ch mkdir /store/group/phys_exotica/dijet/dazsle/TrackObservablesStudy/${ana}`

    # remove existing output files from EOS directory 
    for tfhandle in ${filehandles[*]} ; do
    for tfile in $(xrdfs root://eoscms.cern.ch ls -u /store/group/phys_exotica/dijet/dazsle/TrackObservablesStudy/${ana} | grep ${tfhandle}); do
        echo ${tfile}
        filehandle=(${tfile//\// })
        roothandle=${filehandle[${#filehandle[@]}-1]}
	eval `eos root://cmseos.fnal.gov rm /store/user/cvernier/${ana}/${roothandle}`
        #eval `eos root://eoscms.cern.ch rm /store/group/phys_exotica/dijet/dazsle/TrackObservablesStudy/${ana}/${roothandle}`
    done
        # remove output log files
        prefix="processed-"
        #rm ./${ana}/*${tfhandle#$prefix}*
    done


    # submit condor jobs
    python SubmitCondor.py --indir /store/user/cvernier/fromPhil/WtoQQ/ \
        --outdir ./${ana}/ \
        --maxEvents 50000 \
        --evPerJob 50000 \
        --outdirname ${ana} \
        --cfg ${anacfgs[i]}

    let "i += 1"
done

;;

MERGE )
############################# MERGE #################################
echo "Merging outputs"
cd /uscms_data/d3/${USER}/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd /afs/cern.ch/work/c/cvernier/TrackObservables/trackObservables/processing

for ana in ${anasubs[*]} ; do
for fileh in ${filehandles[*]} ; do
    echo " "
    echo " "
    echo " "
    echo "Working on ${fileh}-${ana}.root"
    echo " "
    hadd -k ${fileh}-${ana}.root  \
        `xrdfs root://eoscms.cern.ch ls -u /store/group/phys_exotica/dijet/dazsle/TrackObservablesStudy/${ana} | grep root | grep ${fileh}`
    mv ${fileh}-${ana}.root ../../samples/
done
done


;;


esac
