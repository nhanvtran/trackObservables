eosdirname=trainings_07152016_all
logdirname=trainings_07152016_all

sigList=(tt WW ZZ)
bkgList=(tt WW ZZ)
treeList=(t_tragam t_tracks t_allpar)

# FULL LIST OF VARIABLES:
# njets 
# j_ptfrac
# j_pt
# j_eta
# j_mass
# j_tau1_b1
# j_tau2_b1
# j_tau3_b1
# j_tau1_b2
# j_tau2_b2
# j_tau3_b2
# j_tau21_b1
# j_tau21_b2
# j_tau21_b1
# j_tau21_b2
# j_tau32_b1
# j_tau32_b2
# j_zlogz
# j_c1_b0
# j_c1_b1
# j_c1_b2
# j_c2_b1
# j_c2_b2
# j_d2_b1
# j_d2_b2
# j_mass_trim
# j_mass_mmdt
# j_mass_prun
# j_mass_sdb2
# j_mass_sdm1
# j_multiplicity

varList=(j_mass_trim j_mass_mmdt j_mass_prun j_mass_sdb2 j_mass_sdm1 j_c2_b1 j_c2_b2 j_d2_b1 j_d2_b2 j_tau21_b1 j_tau21_b2 j_c1_b0 j_c1_b1 j_c1_b2 j_tau1_b1 j_tau1_b2 j_zlogz j_tau32_b1 j_tau32_b2) 


# YOU SHOULD NOT HAVE TO CHANGE ANYTHING BELOW HERE

# setup
SIGS=""
BKGS=""
VARS=""
varsToWeight=(j_mass_trim j_mass_mmdt j_mass_prun j_mass_sdb2 j_mass_sdm1)

# fetch signals
for sig in ${sigList[*]}
do
    if [[ "${SIGS}" == "" ]];
    then
        SIGS="${sig}"
    else
        SIGS="${SIGS},${sig}"
    fi
done

# fetch backgrounds
for bkg in ${bkgList[*]}
do
    if [[ "${BKGS}" == "" ]];
    then
        BKGS="${bkg}"
    else
        BKGS="${BKGS},${bkg}"
    fi
done

# fetch variables to train on
for var in ${varList[*]}
do
    if [[ "${VARS}" == "" ]];
    then
        if [[ "${varsToWeight[@]}" =~ "${var}" ]]; then
            VARS="${var}[0]/j_ptfrac[0]"
        else
            VARS="${var}[0]"
        fi
    else
        if [[ "${varsToWeight[@]}" =~ "${var}" ]]; then
            VARS="${VARS};${var}[0]/j_ptfrac[0]"
        else
            VARS="${VARS};${var}[0]"
        fi
    fi
done

echo "Launching trainings with the following settings:"
echo ""
echo " - SIGNALS: ${SIGS}"
echo " - BACKGROUNDS: ${BKGS}"
echo " - VARIABLES: ${VARS}"
echo ""
echo " - CONDOR OUTPUT DIRECTORY: ./${logdirname}" 
echo " - EOS    OUTPUT DIRECTORY: /eos/uscms/store/user/ecoleman/SubROC/${eosdirname}"
echo ""
echo ""

# make important directories
mkdir /eos/uscms/store/user/ecoleman/SubROC/${eosdirname}
mkdir ${logdirname}

# launch jobs
for tree in ${treeList[*]}
do
    python launchJobs.py -b --doTraining --userOverride ecoleman --tmpDir ${logdirname} \
        --eosDest SubROC/${eosdirname} --treeName ${tree} --sigs ${SIGS} --bkgs ${BKGS} \
        --vars ${VARS}
    
    python launchJobs.py -b --doTraining --userOverride ecoleman --prefix processed-pythia82-fcc100 \
        --ptfix pt5 --tmpDir ${logdirname} --eosDest SubROC/${eosdirname} --treeName ${tree}        \
        --sigs ${SIGS} --bkgs ${BKGS} --vars ${VARS}
done
