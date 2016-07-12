#trackObservables


##installation

Trainings need to be run on cmslpc, and you'll need an EOS area to store condor output. Here is the advised setup:

```
$ cd ~/nobackup/
$ mkdir Substructure-ROC
$ cd Substructure-ROC
$ cmsrel CMSSW_7_2_0
$ git clone https://github.com/nhanvtran/trackObservables.git
$ mkdir /eos/uscms/store/user/${USER}/SubROC/training/weights
```

###analysis

Main body of the BDT training. Scripts to configure and launch TMVA.

You can run trainings via **condor** using `launchJobs.py`, for example:

```
python analysis/launchJobs.py -b --doTraining --makeROCs --treeName t_allpar
```

Input files are assumed to have the form `<prefix>-<type>-<ptfix>-<postfix>.root`. Unless specified otherwise, the defaults for these tags are:

 - prefix : `processed-pythia82-lhc13`, change with `--prefix`
 - type : no default, but can be e.g. `qq`, `tt`, `WW`, `ZZ`, `gg`
 - ptfix : `pt1`, change with `--ptfix`
 - postfix : `50k`, change with `--postfix`

Individual trainings can be run **locally** via `analysis.py`. The following will produce one of the BDTs produced in the condor implementation:

```
nohup python trackObservables/analysis/analysis.py -b \
--sampleDir /uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/prod-Jun14  \
--weightDir ./weights  \
--plotDir   ./plots    \
--sigTag    WW-pt1     \
--bkgTag    ZZ-pt1     \
--inputs    "j_tau21_b1[0];j_tau21_b2[0];j_c1_b0[0];j_c1_b1[0];j_c1_b2[0];j_c2_b1[0];j_c2_b2[0];j_d2_b1[0];j_d2_b2[0];j_mass_trim[0]*j_ptfrac[0];j_mass_mmdt[0]*j_ptfrac[0];j_mass_prun[0]*j_ptfrac[0];j_mass_sdb2[0]*j_ptfrac[0];j_mass_sdm1[0]*j_ptfrac[0];j_mass[0]*j_ptfrac[0]" \
--treeName t_allpar    \
--doTraining > ./output/WW-ZZ-pt1-tracks.txt
```

To reiterate, `--prefix` and `--postfix` will specify the input root files to pull from `--sampleDir`.


To plot ROCs from trainings:

```
$ root -l -b
[0] .L analysis/plotROC.C
[1] plotROC(input1, label1, input2, label2, ..., inputN, labelN, outputName) 
[2] .q
```

To get separation information from BDT trainings (condor only), run `python getSeparationTXT.py` from the directory where you launched the condor jobs. This will output `.txt` files with the desired tables. 

###plotting

Utilities to get formatted control plots from ntuples. Samples available at `/uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/prod-Jun14 `

###processing

Scripts to generate ntuples.



