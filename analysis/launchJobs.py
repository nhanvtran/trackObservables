#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=True, help='no X11 windows')
parser.add_option('--doTraining', action='store_true', dest='doTraining', default=True, help='training or not')
parser.add_option('--cleaning', action='store_true', dest='cleaning', default=True, help='training or not')
parser.add_option('--makeROCs', action='store_true', dest='makeROCs', default=False, help='training or not')
parser.add_option('--userOverride', action='store', type="string", dest='userOverride', default="")
parser.add_option('-i','--interactive', action='store_true', dest='interactive', default=False, help='training or not')
parser.add_option('--treeName',action="store",type="string",dest="treeName",default="t_tragam")
parser.add_option('--sigs',action="store",type="string",dest="sigs",default="tt,WW,ZZ,qq,gg")
parser.add_option('--bkgs',action="store",type="string",dest="bkgs",default="tt,WW,ZZ,qq,gg")
parser.add_option('--prefix',action="store",type="string",dest="prefix",default="processed-pythia82-lhc13")
parser.add_option('--postfix',action="store",type="string",dest="postfix",default="50k")
parser.add_option('--ptfix',action="store",type="string",dest="ptfix",default="pt1")
parser.add_option('--tmpDir',action="store",type="string",dest="tmpDir",default="tmp")
parser.add_option('--eosDest',action="store",type="string",dest="eosDest",default="SubROC/training")
parser.add_option('--vars',action="store",type="string",dest="vars",default="")
parser.add_option('--QGvars',action="store",type="string",dest="QGvars",default="")

parser.add_option('--sampleDir',action="store",type="string",dest="sampleDir",default="/uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/prod-Jun14")

(options, args) = parser.parse_args()


user = os.environ['USER'] if options.userOverride == "" else options.userOverride;

sampleDir = options.sampleDir
#/uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/testdat/';
# weightDir = '/eos/uscms/store/user/ntran/SUSY/theory_JPM/training/weights';
# plotDir   = '/eos/uscms/store/user/ntran/SUSY/theory_JPM/training/plots';
weightDir = './weights/';
eosweightdir = '/store/user/%s/%s/weights/'%(user,options.eosDest);
if options.makeROCs and options.interactive:
    weightDir = './weights';
plotDir   = './plots';

def condorize(command,tag):

    print "--------------------"
    print "Processing sample..."

    prefix = "o";
    prefixR = "R";
    if options.makeROCs: prefix = "op";
    if options.makeROCs: prefixR = "opR";

    # mk eos directory for final ouput
    if not os.path.exists("/eos/uscms/store/user/%s/%s/"%(user,options.eosDest)):
        os.makedirs("/eos/uscms/store/user/%s/%s"%(user,options.eosDest));

    startdir = os.getcwd();
    #change to a tmp dir
    if not os.path.exists(options.tmpDir): os.makedirs(options.tmpDir);
    os.chdir(options.tmpDir);
    curdir = os.getcwd();

    f1n = "tmp_%s.sh" %(tag);
    f1=open(f1n, 'w')
    f1.write("#!/bin/sh \n");

    # setup environment
    f1.write("source /cvmfs/cms.cern.ch/cmsset_default.sh \n");
    f1.write("set SCRAM_ARCH=slc6_amd64_gcc481\n")
    f1.write("cd /uscms_data/d3/%s/Substructure-ROC/CMSSW_7_2_0/src/ \n"%(user))
    f1.write("eval `scramv1 runtime -sh`\n")
    # copy over all stuff and run
    f1.write("cd - \n");
    if options.doTraining: f1.write("mkdir weights \n")
    if options.makeROCs:
        f1.write("mkdir plots \n")
        f1.write("for file in /eos/uscms/%s/o_%s.tar.gz; \n" % (eosweightdir,tag));
        f1.write("do \n");
        f1.write("	filebase=${file##*/} \n");
        f1.write("  echo $filebase \n");
        f1.write("	xrdcp root://cmseos.fnal.gov/%s/$filebase $filebase \n" % (eosweightdir) );
        f1.write("	tar -xvzf $filebase \n");
        f1.write("done \n");

    f1.write("cp %s/*.py . \n" % (startdir));
    f1.write("cp %s/tdrstyle.C . \n" % (startdir));

    f1.write("ls \n");
    f1.write(command+" \n")
    if options.doTraining:
        f1.write("tar -cvzf %s_%s.tar.gz weights \n" % (prefix,tag))
        f1.write("tar -cvzf %s_%s.tar.gz *.root \n" % (prefixR,tag))
        f1.write("xrdcp -f %s_%s.tar.gz root://cmseos.fnal.gov//store/user/%s/%s/%s_%s.tar.gz \n" % (prefix,tag,user,options.eosDest,prefix,tag));
        f1.write("xrdcp -f %s_%s.tar.gz root://cmseos.fnal.gov//store/user/%s/%s/%s_%s.tar.gz \n" % (prefixR,tag,user,options.eosDest,prefixR,tag));
    if options.makeROCs:
        f1.write("tar -cvzf %s_%s.tar.gz plots \n" % (prefix,tag))
        f1.write("xrdcp -f %s_%s.tar.gz root://cmseos.fnal.gov//store/user/%s/%s/%s_%s.tar.gz \n" % (prefix,tag,user,options.eosDest,prefix,tag));
    f1.write("rm *.py *.pyc *.tar.gz *.C *.root \n");
    f1.close()

    f2n = "tmp_%s.condor" % (tag);
    outtag = "out_%s_$(Cluster)" % (tag)
    f2=open(f2n, 'w')
    f2.write("universe = vanilla \n");
    f2.write("Executable = %s \n" % (f1n) );
    f2.write("Requirements = Memory >= 199 &&OpSys == \"LINUX\"&& (Arch != \"DUMMY\" )&& Disk > 1000000 \n");
    f2.write("Should_Transfer_Files = YES \n");
    f2.write("WhenToTransferOutput  = ON_EXIT_OR_EVICT \n");
    f2.write("Output = "+outtag+".stdout \n");
    f2.write("Error = "+outtag+".stderr \n");
    f2.write("Log = "+outtag+".log \n");
    f2.write("Notification    = Error \n");
    f2.write("x509userproxy = $ENV(X509_USER_PROXY) \n")

    f2.write("Queue 1 \n");
    f2.close();

    os.system("condor_submit %s" % (f2n));

    os.chdir("../.");

##################################################################################################
##################################################################################################
##################################################################################################

if __name__ == '__main__':

    ptfix = options.ptfix;

    signalsname = options.sigs.split(',');
    signals =[];
    #masses = ['500','1000','1500','2000'];
    masses = ['1']#,'1000','1500'];
    # masses = ['500','1500'];
    for m in signalsname:
        signals.append( m+"-"+ptfix );

    # vs.
    backgroundsname = options.bkgs.split(',');
    backgrounds = []
    for m in  backgroundsname:
        backgrounds.append( m+"-"+ptfix );

    # observables.append( "mEff" );
    observables = [];
    observablesQG = [];
    if not options.vars=="" :
        observables=options.vars.split(';')
        if options.QGvars=="":
            observablesQG=options.vars.split(';')
        else :
            observablesQG=options.QGvars.split(';')
    else :
        if options.treeName=="t_allpar":
            # mass only
            #observables=["j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
            #observablesQG=["j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
            # shapes only
            observables=["j_tau32_b1[0]","j_tau32_b2[0]","j_tau21_b1[0]","j_tau21_b2[0]","j_c2_b2[0]","j_c2_b1[0]"]
            observablesQG=["j_zlogz[0]","j_tau1_b1[0]","j_tau1_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_multiplicity"]
            # original
            #observables=["j_tau21_b1[0]","j_tau21_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_c2_b1[0]","j_c2_b2[0]","j_d2_b1[0]","j_d2_b2[0]","j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
            #observablesQG=["j_zlogz[0]","j_tau21_b1[0]","j_tau21_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_c2_b1[0]","j_c2_b2[0]","j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
            # tau32
           #observables=["j_tau32_b1[0]","j_tau32_b2[0]","j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
           #observablesQG=["j_tau23_b1[0]","j_tau23_b2[0]","j_mass_trim[0]","j_mass_mmdt[0]","j_mass_prun[0]","j_mass_sdb2[0]","j_mass_sdm1[0]","j_mass[0]"]
        else :
            # mass only
            #observables=["j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]
            #observablesQG=["j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]
            # shapes only
            observables=["j_tau32_b1[0]","j_tau32_b2[0]","j_tau21_b1[0]","j_tau21_b2[0]","j_c2_b2[0]","j_c2_b1[0]"]
            observablesQG=["j_zlogz[0]","j_tau1_b1[0]","j_tau1_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_multiplicity"]
            # original
            #observables=["j_tau21_b1[0]","j_tau21_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_c2_b1[0]","j_c2_b2[0]","j_d2_b1[0]","j_d2_b2[0]","j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]
            #observablesQG=["j_zlogz[0]","j_tau21_b1[0]","j_tau21_b2[0]","j_c1_b0[0]","j_c1_b1[0]","j_c1_b2[0]","j_c2_b1[0]","j_c2_b2[0]","j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]
            # tau32
            # observables=["j_tau32_b1[0]","j_tau32_b2[0]","j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]
             #observablesQG=["j_zlogz[0]","j_tau23_b1[0]","j_tau23_b2[0]","j_mass_trim[0]*j_ptfrac[0]","j_mass_mmdt[0]*j_ptfrac[0]","j_mass_prun[0]*j_ptfrac[0]","j_mass_sdb2[0]*j_ptfrac[0]","j_mass_sdm1[0]*j_ptfrac[0]","j_mass[0]*j_ptfrac[0]"]



    ###########
    ## training
    jobctr = 0;
    Observables = []
    if options.doTraining:
        for sig in signals:
            for bkg in backgrounds:
                obsList = '';
                if sig == bkg: continue
                if (sig == "qq-%s"%(ptfix) or sig == "gg-%s"%(ptfix)) and (bkg == "qq-%s"%(ptfix) or bkg == "gg-%s"%(ptfix)):
                    Observables = observablesQG
                else : Observables = observables
                for iObs in range(len(Observables)):
                    if iObs != len(Observables)-1:
                        tmp = Observables[iObs]
                        #tmp2= tmp.replace('[','')
                        #tmp3 = tmp2.replace(']','')
                        obsList += tmp + ";";
                    else:
                        tmp = Observables[iObs]
                        #tmp2= tmp.replace('[','')
                        #tmp3 = tmp2.replace(']','')
                        obsList += tmp;

                command = "python analysis.py -b ";
                command += " --sampleDir " + sampleDir;
                command += " --weightDir " + weightDir;
                command += " --plotDir "   + plotDir;
                command += " --sigTag "    + sig;
                command += " --bkgTag "    + bkg;
                command += " --treeName "   + options.treeName;
                command += " --prefix "    + options.prefix;
                command += " --postfix "   + options.postfix;
                command += " --userOverride "   + options.userOverride;
                command += " --inputs "         + '"'+obsList+'"';
                command += " --doTraining";

                print command;
                label = sig + "_" + bkg+ "_"+ options.treeName;
                #tag = label.translate(None, ';,');
                labelbase = sig + "_" + bkg+ "_"+ options.treeName;
                tag = label.replace(',','_');
                tagbase = labelbase.replace(',','_')

                filestring = "%s/weights_bdtg_%s/TMVAClassification_MVA_bdtg_%s_BDTG.weights.xml" % (weightDir,tagbase,tag);
                print("here")
                if options.interactive:
                    print "doing, ", command
                    os.system(command);
                    print("here")
                else:
                    if not os.path.isfile(filestring) and options.cleaning:
                        jobctr+=1;
                        condorize(command,tag);
                        print("here1")
                    if not options.cleaning:
                        jobctr+=1;
                        condorize(command,tag);
                        print("here2")
                    time.sleep(0.1) # delays for 5 seconds

    ## make roc
    if options.makeROCs:
        for sig in signals:
            for bkg in backgrounds:
                obsList = '';
                if sig == bkg: continue
                if (sig == "qq-%s"%(ptfix) or sig == "gg-%s"%(ptfix)) and (bkg == "qq-%s"%(ptfix) or bkg == "gg-%s"%(ptfix)) :
                    Observables = observablesQG
                else : Observables = observables
                for iObs in range(len(Observables)):
                    if iObs != len(Observables)-1: obsList += Observables[iObs] + ";";
                    else: obsList += Observables[iObs];

                command = "python analysis.py -b ";
                command += " --sampleDir " + sampleDir;
                command += " --weightDir " + weightDir;
                command += " --plotDir "   + plotDir;
                command += " --sigTag "    + sig;
                command += " --bkgTag "    + bkg;
                command += " --treeName "   + options.treeName;
                command += " --prefix "    + options.prefix;
                command += " --postfix "   + options.postfix;
                command += " --userOverride "   + options.userOverride;
                command += " --inputs "    + '"'+obsList+'"';
                command += " --makeROCs";

                print command;
                label = sig + "_" + bkg + "_ROCs";
                #tag = label.translate(None, ';,');
                labelbase = sig + "_" + bkg;
                tag = label.replace(',','_');
                tagbase = labelbase.replace(',','_')

                filestring = "%s/plots_bdtg_%s/RocSummary_%s.txt" % (plotDir,tagbase,tagbase);
                #print filestring
                # print tag;
                # condorize(command,tag);
                if options.interactive:
                    print "interactively: ", command, tagbase
                    os.system(command);
                else:
                    if not os.path.isfile(filestring) and options.cleaning:
                        jobctr+=1;
                        condorize(command,tagbase);
                    if not options.cleaning:
                        jobctr+=1;
                        condorize(command,tagbase);
                    time.sleep(0.1) # delays for 5 seconds

    print "total jobs = ", jobctr;
