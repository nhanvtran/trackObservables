#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from observableContainer import *
from TMVAhelper import *
from utilities import *
from subprocess import call

ROOT.gROOT.ProcessLine(".L tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadTopMargin(0.06);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

parser.add_option('--doTraining', action='store_true', dest='doTraining', default=True, help='perform training')
parser.add_option('--makeROCs', action='store_true', dest='makeROCs', default=False, help='produce ROC curves')

parser.add_option('--inputs',action="store",type="string",dest="inputs",default="HT")
parser.add_option('--sampleDir',action="store",type="string",dest="sampleDir",default="/uscms_data/d2/ntran/physics/FCC/trackObservablesStudy/trackObservables/processing/prod-Jun14")
parser.add_option('--sigTag',action="store",type="string",dest="sigTag",default="qq-pt5")
parser.add_option('--bkgTag',action="store",type="string",dest="bkgTag",default="gg-pt5")
parser.add_option('--weightDir',action="store",type="string",dest="weightDir",default="/store/user/${USER}/SubROC/training/weights")
parser.add_option('--userOverride',action="store",type="string",dest="userOverride",default="")
parser.add_option('--prefix',action="store",type="string",dest="prefix",default="processed-pythia82-lhc13")
parser.add_option('--postfix',action="store",type="string",dest="postfix",default="50k")
parser.add_option('--plotDir',action="store",type="string",dest="plotDir",default="/store/user/${USER}/SubROC/training/plots")
parser.add_option('--treeName',action="store",type="string",dest="treeName",default="t_tragam")

# TODO: needed? can be replaced by --postfix
parser.add_option('--trainingSample', action='store_true', dest='trainingSample', default=False, help='training or not')

(options, args) = parser.parse_args()

############################################################

# observableTraining takes in 2 root files, list of observables, spectator observables ... launches a CONDOR job
# TMVAhelper.py is used by observableTraining
# analysis.py defines the list of trainings, the signal and background process

########################################################################################################################
########################################################################################################################

if __name__ == '__main__':

    user = os.environ['USER'] if options.userOverride == "" else options.userOverride;

    prefix = options.prefix;
    postfix = options.postfix;
	#if (options.sigTag == 'ZZ-pt5' or options.bkgTag == 'ZZ-pt5'): postfix = '10k';
	#if options.trainingSample: postfix = 'train'

    weightloc = options.weightDir.replace("${USER}", user);
    odir = options.plotDir.replace("${USER}", user);
    sampleDir = options.sampleDir.replace("${USER}", user);
    sigName = options.sigTag;
    bkgName = options.bkgTag;
    treeName = options.treeName;
    sigFN = sampleDir+"/"+prefix+"-"+sigName+"-"+postfix+".root";
    bkgFN = sampleDir+"/"+prefix+"-"+bkgName+"-"+postfix+".root";

    variables = [];
    for ivars in options.inputs.split(';'):
        variables.append( ivars.split(',') );
    print sigFN, bkgFN, variables

    f_sig = ROOT.TFile(sigFN);
    f_bkg = ROOT.TFile(bkgFN);
    t_sig = f_sig.Get(treeName);
    t_bkg = f_bkg.Get(treeName);
    print(t_sig.GetEntries())
    print(t_bkg.GetEntries())

    cuts = [];
    cuts.append( ["j_passCut",1,1] );
#	cuts.append( ["HT",500,99999] );
#	cuts.append( ["MHT",200,99999] );
#	cuts.append( ["NJets",1,99] );

#                   NAME        MIN   MAX
    MVAvariables = { "njets":    (0 , 9),
            "j_ptfrac[0]":       (0 , 1),
            "j_mass[0]":         (0 , 500),
            "j_tau1_b1[0]":      (0 , 300),
            "j_tau2_b1[0]":      (0 , 300),
            "j_tau1_b2[0]":      (0 , 300),
            "j_tau2_b2[0]":      (0 , 300),
            "j_tau21_b1[0]":     (0 , 1.5),
            "j_tau21_b2[0]":     (-5, 1.5),
            "j_zlogz[0]":        (-5, 5),
            "j_c1_b0[0]":        (-5, 5),
            "j_c1_b1[0]":        (-5, 5),
            "j_c1_b2[0]":        (-5, 5),
            "j_c2_b1[0]":        (-5, 5),
            "j_c2_b2[0]":        (0 , 5),
            "j_d2_b1[0]":        (0 , 100),
            "j_d2_b2[0]":        (0 , 100),
            "j_mass_trim[0]":    (0 , 500),
            "j_mass_mmdt[0]":    (0 , 500),
            "j_mass_prun[0]":    (0 , 500),
            "j_mass_sdb2[0]":    (0 , 500),
            "j_mass_sdm1[0]":    (0 , 500),
            "j_multiplicity[0]": (0 , 120) }

    allVariables = [(key) for key in MVAvariables]
    varlo = [(MVAvariables[key][0]) for key in MVAvariables]
    varhi = [(MVAvariables[key][1]) for key in MVAvariables]

    h_variables = [];
    for i in range(len(allVariables)):
        h_variables.append( [ROOT.TH1F("hs_"+allVariables[i],"au; "+allVariables[i]+";",20,varlo[i],varhi[i]), ROOT.TH1F("hb_"+allVariables[i],"au; "+allVariables[i]+";",20,varlo[i],varhi[i]) ] );

    # variables = [];
    # variables.append( ["HT"] );
    # # variables.append( ["NJets"] );
    # # variables.append( ["MHT"] );
    # # variables.append( ["HT","MHT"] );
    # # variables.append( ["HT","NJets"] );
    # # variables.append( ["NJets","MHT"] );
    # # variables.append( ["HT","NJets","MHT"] );
    # spectators = [];

    observableSets = [];
    h_bdts = [];
    tagbase = "%s_%s_%s" % (sigName,bkgName,treeName);
    labelbase = "bdtg_%s_%s_%s" % (sigName,bkgName,treeName);

    # make ROC dir
    curodir = odir + "/" + "plots_" + labelbase;
    if options.makeROCs:
        if not os.path.exists(curodir): os.makedirs(curodir);

    #for i in range(len(variables)):
    #	for vnames in variables[i]: label += "_" + vnames;

    # make weights dir
    curweightloc = weightloc+"/" + "weights_" + labelbase;
    os.makedirs(curweightloc)
    if options.doTraining:
        if not os.path.exists(curweightloc): os.makedirs(curweightloc);
            # currentFilesInDir = os.listdir(curweightloc);
            # if label in currentFilesInDir: continue;

    # perform training, if requested
    observableSets.append( observableContainer(sigFN,bkgFN,variables,cuts,labelbase,options.treeName,curweightloc) );
    if options.doTraining: observableSets[0].doTraining();
    observableSets[0].readMVA("BDTG");
    h_bdts.append( [ROOT.TH1F("hsbdt_"+labelbase,"au; "+labelbase+";",100000,-1,1), ROOT.TH1F("hbbdt_"+labelbase,"au; "+labelbase+";",100000,-1,1) ] );

    ##### ##### ##### ##### #####
    # FILL TREES
    if options.makeROCs:

        trees = [t_sig,t_bkg];
        names = [sigName,bkgName];
        for it in range(len(trees)):
            nent = trees[it].GetEntriesFast()
            for i in range(nent):

                # if i > 1000: break;

                if(i % (1 * nent/100) == 0):
                    sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
                    sys.stdout.flush();


                trees[it].GetEntry(i);

                # print "lheWeight = ", trees[it].lheWeight

                passesCuts = True;
                for cut in cuts:
                    if getattr(trees[it],cut[0]) < cut[1] or getattr(trees[it],cut[0]) > cut[2]: passesCuts = False;
                if passesCuts == False: continue;

                # fill the variable histograms
                for ivar in range(len(allVariables)):
                    # print getattr(trees[it],allVariables[ivar]), getattr(trees[it],"lheWeight") ;
                    h_variables[ivar][it].Fill( getattr(trees[it],allVariables[ivar]), getattr(trees[it]) );
                # fill the bdt histograms
                #for ibdt in range(len(variables)):
                tmplist = [];
                #	for var in variables[ibdt]: tmplist.append( getattr(trees[it],var) );
                h_bdts[it].Fill( observableContainer.evaluateMVA(tmplist,"BDTG"), getattr(trees[it]) );

            print "\n";

    ##### ##### ##### ##### #####
    # PLOT
        print "make raw variable plots..."
        for ivar in range(len(allVariables)):
            makeCanvas(h_variables[ivar],names,allVariables[ivar]+'_'+tagbase,curodir,True);
        print "make bdt distribution plots..."
        #for ibdt in range(len(variables)):
        label = "bdt_";
        #	for vnames in variables[ibdt]: label += vnames + "_";
        #	print label;
        makeCanvas(h_bdts,names,label+'_'+tagbase,curodir,True,True);

    ##### ##### ##### ##### #####
    # make ROCs
        print "make rocs..."
        rocs = [];
        leg = ROOT.TLegend( 0.2, 0.6, 0.5, 0.9 );
        leg.SetBorderSize( 0 );
        leg.SetFillStyle( 0 );
        leg.SetTextSize( 0.03 );

        labels2 = [];
        #for ibdt in range(len(variables)):
        print "roc #";
        rocs.append( makeROCFromHisto(h_bdts) );
        rocs.SetLineColor(1);
        rocs.SetLineWidth(2)
        label = "";
        label2 = "";
        #for ivar in range(len(variables[ibdt])):
        #		label += str(variables[ibdt][ivar]);
    #			label2 += str(variables[ibdt][ivar]);
    #			if ivar < len(variables[ibdt])-1:
    #				label += "+";
    #				label2 += "_"

        leg.AddEntry( rocs, label, "l" );
        labels2.append(label2);

        bkgrej = [];

        canroc = ROOT.TCanvas("canroc","canroc",1200,1000);
        hrl1 = canroc.DrawFrame(5e-2,1e-6,1.0,1.0);
        # hrl1.GetXaxis().SetTitle("#varepsilon_{sig} ("+signame+")");
        # hrl1.GetYaxis().SetTitle("#varepsilon_{bkg} ("+bkgname+")");
        hrl1.GetXaxis().SetTitle("signal efficiency");
        hrl1.GetYaxis().SetTitle("background efficiency");
        for i in range(len(rocs)):
            rocs[i].Draw("l");
            tmpbkgrej = [];
            tmpbkgrej.append( labels2[i] );
            if rocs[i].Eval(0.1) > 0: tmpbkgrej.append( 1./rocs[i].Eval(0.1) );
            else: tmpbkgrej.append( -1 );
            if rocs[i].Eval(0.25) > 0: tmpbkgrej.append( 1./rocs[i].Eval(0.25) );
            else: tmpbkgrej.append( -1 );

            bkgrej.append( tmpbkgrej );

        leg.Draw();
        canroc.SaveAs(curodir+"/Rocs"+'_'+tagbase+".root");
        canroc.SaveAs(curodir+"/Rocs"+'_'+tagbase+".eps");
        canroc.SaveAs(curodir+"/Rocs"+'_'+tagbase+".png");
        canroc.SaveAs(curodir+"/Rocs"+'_'+tagbase+".pdf");
        ROOT.gPad.SetLogy();
        ROOT.gPad.SetLogx();
        canroc.SaveAs(curodir+"/Rocs_log"+'_'+tagbase+".root");
        canroc.SaveAs(curodir+"/Rocs_log"+'_'+tagbase+".eps");
        canroc.SaveAs(curodir+"/Rocs_log"+'_'+tagbase+".png");
        canroc.SaveAs(curodir+"/Rocs_log"+'_'+tagbase+".pdf");

        fout = open(curodir+"/RocSummary_"+tagbase+".txt",'w');
        for line in bkgrej:
            ostring = '';
            for a in range(len(line)): ostring += str(line[a]) + ' ';
            ostring += '\n';
            fout.write(ostring)
        fout.close();

