#!/bin/bash

if [[ "$1" == "" ]]; then
    echo "*****************************************"
    echo "*                                       *"
    echo "* Track Observables Study: Steer Script *"
    echo "*                                       *"
    echo "*****************************************"
    echo "* MAKE                                  *"
    echo "* JOBS                                  *"
    echo "* MERGE                                 *"
    echo "* PLOT                                  *"
    echo "* WWW                                   *"
    echo "*****************************************"
    echo ""
    
    exit 1
fi

procdir=$PWD

anasubs=(nores_nogran res_nogran res_granH res_granHE resx100_granHE)
anacfgs=(z r rh rhe rhex)

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



####################################################################
#             Nothing should change below this line                #
####################################################################


case $1 in

MAKE )
############################# MAKE #################################
echo "Making executables"

make
eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/TrackObservablesStudy/anaSubstructure`
xrdcp anaSubstructure root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/

;;

JOBS)
############################# JOBS #################################
echo "Running jobs"

i=0
for ana in ${anasubs[*]} ; do
    # create output directory, if it does not exist
    eval `eos root://cmseos.fnal.gov mkdir /store/user/${USER}/${ana}`

    # remove existing output files from EOS directory 
    for tfile in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/${ana}); do
        filehandle=(${tfile//\// })
        roothandle=${filehandle[${#filehandle[@]}-1]}
        eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/${ana}/${roothandle}`
    done

    # remove output log files
    rm ./${ana}/*

    # submit condor jobs
    python SubmitCondor.py --indir /store/user/${USER}/TrackObservablesStudy/fromMarat/ \
        --outdir ./${ana}/ \
        --maxEvents 50000 \
        --evPerJob 10000 \
        --anaSubLoc root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/anaSubstructure \
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

cd /uscms_data/d3/${USER}/TrackObservablesStudy/trackObservables/processing/

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

;;


PLOT )
############################## PLOT #################################

cd /uscms_data/d3/ecoleman/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd ${procdir}/../plotting/

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


;;

WWW )
############################### WWW #################################
echo "Sending plots to your www directory on LXPLUS"

cd ${procdir}/../plotting/

if [[ "$2" == "" ]]; then
    echo ""
    echo "WWW requires a second argument for the target directory in your www area."
    echo " - format: ~/www/TrackObservablesStudy/SmearPlots/\$2/"
    exit 1;
fi

for ana in ${anasubs[*]} ; do
    scp plots/${ana}/* ${USER}@lxplus.cern.ch:~/www/TrackObservablesStudy/SmearPlots/$2/${ana}/
done

;;

esac
