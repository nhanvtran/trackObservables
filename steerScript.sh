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
    echo "* PLOT   - plot output                  *"
    echo "* WWW    - send plots to CERNWEB        *"
    echo "* TRAIN  - run BDT trainings            *"
    echo "* BDTS_MOVE - move BDT trainings        *"
    echo "* BDTS_PLOT - plot BDT results          *"
    echo "* BDTS_WWW  - send BDT plots to CERNWEB *"
    echo "*****************************************"
    echo ""
    
    exit 1
fi

cd processing
procdir=$PWD

ntupleLoc=/store/user/ecoleman/TrackObservablesStudy/fromMarat/ 

anasubs=(r1_h022_e0175_t220_nonu
         r0_h0_e0_nonu)
anacfgs=(r0.5:h0.010:e0.005:t220.0
         z)

#anasubs=(r05_h01_e005_t220_nonu
#         r05_h002_e001_t220_nonu
#         r1_h022_e0175_t220_nonu
#         r0_h0_e0_nonu)
#anacfgs=(r0.5:h0.010:e0.005:t220.0
#r0.5:h0.002:e0.001:t220.0 
#r1.0:h0.022:e0.0175:t220.0
#z)

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

trainings=(all
massonly
shapesonly)


####################################################################
#             Nothing needs to change below this line              #
####################################################################
case $1 in

MAKE )
############################# MAKE #################################
echo "Making executables"

make
eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/TrackObservablesStudy/anaSubstructure`
xrdcp -f anaSubstructure root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/
xrdcp -f lhc14-pythia8-4C-minbias-nev100.pu14.gz root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/

shopt -s expand_aliases
alias cachedir='echo "Signature: 8a477f597d28d172789f06886806bc55\n# This file is a cache directory tag.\n# For information about cache directory tags, see:\n#       http://www.brynosaurus.com/cachedir/" > CACHEDIR.TAG'

current=`pwd`
cd $CMSSW_BASE/tmp
cachedir
cd $CMSSW_BASE/src/.git
cachedir
cd $CMSSW_BASE/src/TreeMaker/.git
cachedir
cd $current

tar --exclude-caches-all -zcf ${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION}
xrdcp -f ${CMSSW_VERSION}.tar.gz root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/

;;

JOBS)
############################# JOBS #################################
echo "Running jobs"

i=0
for ana in ${anasubs[*]} ; do
    # create output directory, if it does not exist
    eval `eos root://cmseos.fnal.gov mkdir /store/user/${USER}/${ana}`

    # remove existing output files from EOS directory 
    for tfhandle in ${filehandles[*]} ; do
    for tfile in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/${ana} | grep ${tfhandle}); do
        echo ${tfile}
        filehandle=(${tfile//\// })
        roothandle=${filehandle[${#filehandle[@]}-1]}
        eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/${ana}/${roothandle}`
    done
        # remove output log files
        prefix="processed-"
        rm ./${ana}/*${tfhandle#$prefix}*
    done


    # submit condor jobs
    python SubmitCondor.py --indir ${ntupleLoc} \
        --outdir ./${ana}/ \
        --maxEvents 5000 \
        --evPerJob 500 \
        --anaSubLoc root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/ \
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
    hadd -f -k ${fileh}-${ana}.root  \
        `xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/${ana} | grep root | grep ${fileh}`
    mv ${fileh}-${ana}.root samples/.

done
done
;;

ADD_CUTS )
############################# ADD CUTS #################################
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

    # add in various cuts to samples depending on pt of jet
    if [[ "${fileh}" =~ "pt1" ]] ; then 
        root -l -b -q addCutsToSamples.cxx"(\"../../samples/${fileh}-${ana}.root\",800,1600)"
    elif [[ "${fileh}" =~ "pt5" ]] ; then
        root -l -b -q addCutsToSamples.cxx"(\"../../samples/${fileh}-${ana}.root\",4500,6000)"
    fi

done
done

;;

PLOT )
############################## PLOT #################################

cd /uscms_data/d3/${USER}/CMSSW_8_0_4/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

cd ${procdir}/../plotting/

mkdir plots
for ana in ${anasubs[*]} ; do
    mkdir plots/${ana}
    rm plots/${ana}/*
    
    cmd=""
for fileh in ${filehandles[*]} ; do
    echo " "
    echo " "
    echo " "
    echo "Working on ${fileh}-${ana}.root"
    echo " "

    cmd="${cmd} && python quickPlotter.py -b"
    cmd="${cmd} --basedir ../../samples"
    cmd="${cmd} --ana ${ana}"
    cmd="${cmd} -o ./plots/${ana}/ > ${fileh}-${ana}.out"
done

    nohup sh -c "eval ${cmd}" &
done

;;

WWW )
############################### WWW #################################
echo "Sending plots to your www directory on LXPLUS: ~/www/TrackObservablesStudy/$2/"

cd ${procdir}/../plotting/plots

if [[ "$2" == "" ]]; then
    echo ""
    echo "WWW requires a second argument for the target directory in your www area."
    echo " - format: ~/www/TrackObservablesStudy/\$2/"
    exit 1;
fi


rm plots.tar.gz
find ./*/*.{png,pdf,C} > filenames.txt 
tar -czvf plots.tar.gz -T filenames.txt 
scp plots.tar.gz \
    ${USER}@lxplus.cern.ch:~/www/TrackObservablesStudy/$2/

;;

TRAIN )
############################# TRAIN #################################
echo "Running trainings"

# signal process groups to train (csv)
aloSigList=(WW,ZZ
tt
qq
WW,ZZ
WW)
#aloSigList=(WW)

# background process groups to train (csv)
aloBkgList=(gg,qq
gg,qq
gg
tt
ZZ)
#aloBkgList=(qq)

# variables for each group 
massVars="j_mass_trim,j_mass_mmdt,j_mass_prun,j_mass_sdb2,j_mass_sdm1,j_mass"

aloVarList=("j_c2_b1,j_c2_b2,j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2"
"j_tau32_b1,j_tau32_b2"
"j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_multiplicity"
"j_c2_b1,j_c2_b2,j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2,j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_tau32_b1,j_tau32_b2"
"j_c2_b1,j_c2_b2,j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2,j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_tau32_b1,j_tau32_b2")

shapeVarList=("j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2"
"j_tau32_b1,j_tau32_b2"
"j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_multiplicity"
"j_c2_b1,j_c2_b2,j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2,j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_tau32_b1,j_tau32_b2"
"j_c2_b1,j_c2_b2,j_d2_b1,j_d2_b2,j_tau21_b1,j_tau21_b2,j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_tau32_b1,j_tau32_b2")


cd ${procdir}/../analysis/
for ana in ${anasubs[*]} ; do
    for training in ${trainings[*]} ; do
        echo "Removing files for ${ana} ${training}"
        rm ${procdir}/../analysis/trainings_${ana}_${training}/*

        for tfile in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/SubROC/trainings_${ana}_${training}/); do
            echo ${tfile}
            filehandle=(${tfile//\// })
            roothandle=${filehandle[${#filehandle[@]}-1]}
            eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/SubROC/trainings_${ana}_${training}/${roothandle}`
        done
    done

    i=0
    echo "Launching jobs for ${ana}"
    for sig in ${aloSigList[*]}; do
        
        # shapesonly trainings
        if [[ "${trainings[@]}" =~ "shapesonly" ]] ; then
            cmd="sh launchJobs.sh ${ana} ${aloSigList[$i]} ${aloBkgList[$i]} ${shapeVarList[$i]} shapesonly"
            eval "$cmd" &
        fi

        # massonly trainings
        if [[ "${trainings[@]}" =~ "massonly" ]] ; then
            cmd="sh launchJobs.sh ${ana} ${aloSigList[$i]} ${aloBkgList[$i]} ${massVars} massonly"
            eval "$cmd" &
        fi
        
        # all trainings
        cmd="sh launchJobs.sh ${ana} ${aloSigList[$i]} ${aloBkgList[$i]} ${massVars},${aloVarList[$i]} all"
        eval "$cmd"

        sleep 10

        let "i += 1"    
    done
done

;;

BDTS_MOVE )
############################ BDT_MOVE ################################
echo "Moving BDT training files from EOS"
echo ""

for ana in ${anasubs[*]} ; do
for training in ${trainings[*]} ; do
    echo "For ${ana}_${training}:"
    mkdir -p ${procdir}/../output/${ana}/eosrootfiles/${training}
    rm   -rf ${procdir}/../output/${ana}/eosrootfiles/${training}/*
        
    # move existing output files from EOS directory 
    echo " - Moving files..."
    for tfile in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/SubROC/trainings_${ana}_${training}); do
        eval `xrdcp ${tfile} ${procdir}/../output/${ana}/eosrootfiles/${training}` 
    done
    
    # untar files 
    echo " - Untarring files..."
    cd ${procdir}/../output/${ana}/eosrootfiles/${training}/
    for tfile in $(ls -u ${procdir}/../output/${ana}/eosrootfiles/${training}/*.gz); do
        tar -xvf ${tfile}
    done

    rm     ${procdir}/../output/${ana}/eosrootfiles/${training}/*.gz
    rm -rf ${procdir}/../output/${ana}/eosrootfiles/${training}/weights
done
done

;;

BDTS_PLOT )
############################ BDT_PLOT ################################
echo "Plotting ROCs and grabbing separation statistics"
echo ""

cd ${procdir}/../analysis/

for ana in ${anasubs[*]} ; do
for training in ${trainings[*]} ; do
    mkdir -p ${procdir}/../output/${ana}/${training}
    echo " - Plotting for ${ana}_${training}"
    cd ${procdir}/../analysis/

    ### get separation statisics 
    python getSeparationTXT.py \
        --inputs ./trainings_${ana}_${training}/ \
        --output ${procdir}/../output/${ana}/${training}

    ### get ROC background rejection grids 
    python getROCBkgRej.py \
        --inputs ${procdir}/../output/${ana}/eosrootfiles/${training}/ \
        --output ${procdir}/../output/${ana}/${training} \
        --ana ${ana} --training ${training}
done
    # get ROC plots 
    root -l -b -q "getROCs.cc(\"${procdir}/../output/${ana}/eosrootfiles\",\"${procdir}/../output/${ana}\",\"${ana}\")" 

done

### get ROC background rejection grids (Detector Comparison) 
python getDetecCompBkgRej.py \
    --inputs "${procdir}/../output/r*t*/eosrootfiles/allfin/*.root" \
    --output ${procdir}/../output/ \
    --tree allpar --training allfin10 \
    --name DetecComp_bkgRejPwr_at_sig0p5_

root -l -b -q "getROCsForAllSmearings.cc(\"${procdir}/../output/\",\"${procdir}/../output/\")" & 

;;

BDTS_WWW )
############################ BDT_WWW ################################
echo "Sending BDT plots to your www directory on LXPLUS: "

if [[ "$2" == "" ]]; then
    echo ""
    echo "WWW requires a second argument for the target directory in your www area."
    echo " - format: ~/www/TrackObservablesStudy/\$2/"
    exit 1;
fi

rm output/output.tar.gz
cd ${procdir}/../output/
#tar -czvf output.tar.gz  ./Summary* --exclude=*.{root,lhe,gz} --exclude=*eosrootfiles*
#tar -czvf output.tar.gz  ./Summary* ./*cmpsmear* ./ROCEnvelope* --exclude=*.{root,lhe,gz} --exclude=*eosrootfiles*
#tar -czvf output.tar.gz  ./*/*/BackgroundRejection* --exclude=*.{root,lhe,gz} --exclude=*eosrootfiles*
tar -czvf output.tar.gz ./r0_h0_e0/*.{pdf,png,C} ./r0_h0_e0/*/*.{pdf,png,C} --exclude=*.{root,lhe,gz} --exclude=*eosrootfiles*
scp ${procdir}/../output/output.tar.gz \
    ${USER}@lxplus.cern.ch:~/www/TrackObservablesStudy/$2/

;;

SUMMARY )
############################ SUMMARY ################################
echo "Making summary projections"

# hardcoded settings: anaSubs to use
aloAna="r1_h022_e0175_t220_nonu,r05_h01_e005_t220_nonu,r05_h002_e001_t220_nonu"
#aloAna="r0_h0_e0_nonu"

# kinds of plots:
# 1) a plot for every (tree,pt); draw (sig,ana)
# 2) a plot for every (pt,ana) ; draw (sig,tree)
# 3) a plot for every (pt,ana) ; draw (sig,tree)
lines=(
*,pt1,**,*)
# a list of csv of signals to use for each plot
sigs=(W,Z) 
# a list of csv of variables for each set of plots
varLists=(
j_mass_mmdt
)
#j_c1_b0,j_c1_b1,j_c1_b2,j_tau1_b1,j_tau1_b2,j_zlogz,j_multiplicity,j_mass_trim,j_mass_mmdt,j_mass_prun,j_mass_sdb2,j_mass_sdm1,j_mass
# name indicating kind of plot
nameList=(WvZDetecComp)

cd ${procdir}/../plotting/

i=0
cmd="echo 'plotting...'"
for sig  in ${sigs[*]} ; do

    cmd="${cmd} && python plotSampleSummary.py"
    cmd="${cmd} --basedir /uscms_data/d3/ecoleman/TrackObservablesStudy/samples/"
    cmd="${cmd} --sigs ${sig}"
    cmd="${cmd} --logPlots"
    cmd="${cmd} --vars ${varLists[i]}"
    cmd="${cmd} --ana ${aloAna}"
    cmd="${cmd} --lines ${lines[$i]}"
    cmd="${cmd} -n ${nameList[$i]}"
    cmd="${cmd} -o ${procdir}/../output/ "
    
    let "i+=1"

done
    
eval "${cmd}"

;;
ANATEST )
############################ ANATEST ################################
testconfig="z"

echo "Test anasubstructure with config ${testconfig}"

cd processing

xrdcp root://cmseos.fnal.gov:///store/user/ecoleman/TrackObservablesStudy/fromMarat/pythia82-lhc13-WW-pt1-50k-4.lhe .
./anaSubstructure pythia82-lhc13-WW-pt1-50k-4 ./ 200 500 ${testconfig} 1000 

;;
esac
