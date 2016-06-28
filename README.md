#trackObservables


##installation

###analysis

Main body of the BDT training. Scripts to configure and launch TMVA.

You can run trainings via **condor** using `launchJobs.py`, for example:

```
python analysis/launchJobs.py -b --doTraining --makeROCs --treeName t_allpar
```

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


To plot ROCs from trainings:

```
$ root -l
[0] .L analysis/plotROC.C
[1] plotROC(input1, label1, input2, label2, ..., inputN, labelN) 
[2] .q
```

###plotting

Utilities to get formatted control plots from ntuples. Samples available at `/uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/prod-Jun14 `

###processing

Scripts to generate ntuples.



