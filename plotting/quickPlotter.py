#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import string
import ROOT

import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.12);
ROOT.gStyle.SetPadRightMargin(0.25);
ROOT.gStyle.SetPaintTextFormat("1.1f");

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--mixCategories', action='store_true', dest='mixCats', default=False, help='no X11 windows')
parser.add_option('--logPlots', action='store_false', dest='makeLogPlots', default=True, help='make log plots')
parser.add_option('-i', action='store', dest='files', default="", help='file list of inputs ("parton,file,category;...")')
parser.add_option('-o', action='store', dest='outdir',default="./plots/", help='where to store plotting output')

(options, args) = parser.parse_args()

############################################################

# observableTraining takes in 2 root files, list of observables, spectator observables ... launches a CONDOR job
# TMVAhelper.py is used by observableTraining
# analysis.py defines the list of trainings, the signal and background process

########################################################################################################################
########################################################################################################################

files=[('q','../testSamples/processed-pythia82-lhc13-qq-pt1-50k.root','pt1'),
        ('g','../testSamples/processed-pythia82-lhc13-gg-pt1-50k.root','pt1'),
        ('t','../testSamples/processed-pythia82-lhc13-tt-pt1-50k.root','pt1'),
        ('W','../testSamples/processed-pythia82-lhc13-WW-pt1-50k.root','pt1'),
        ('Z','../testSamples/processed-pythia82-lhc13-ZZ-pt1-50k.root','pt1'),
        ('q','../testSamples/processed-pythia82-fcc100-qq-pt5-50k.root','pt5'),
        ('g','../testSamples/processed-pythia82-fcc100-gg-pt5-50k.root','pt5'),
        ('t','../testSamples/processed-pythia82-fcc100-tt-pt5-50k.root','pt5'),
        ('W','../testSamples/processed-pythia82-fcc100-WW-pt5-50k.root','pt5'),
        ('Z','../testSamples/processed-pythia82-fcc100-ZZ-pt5-50k.root','pt5')] if options.files=="" else [(x.split(',')[0],x.split(',')[1],x.split(',')[2]) for x in options.files.split(';')];

# MECHANISM FOR MODIFYING OPTIONS FOR LOTS OF PLOTS
# - keep in mind that the order matters and is not necessarily the one you input
#    (thanks python), so be as specific as possible always
# - can use regexp's for all categories
newPlotLimits= {
#   FORMAT:
#   "p1;;p2;;ptCat1;;ptCat2;;tree;;branch" : (nBins,min,max),
    # TAU 1,2,3 B1
    "[W,Z,q];;[W,Z,q];;pt;;pt;;[tracks,allpar,tragam];;tau[1,2,3]_b1": (150,0,150),
    "[t,g];;[W,t,g,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;tau1_b1": (250,0,250),
    "[t];;[W,t,g,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau1_b1": (250,0,250),
    "[g];;[W,t,g,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau1_b1": (150,0,150),

    "[W,q,Z];;[W,g,t,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau2_b1": (70,0,70),
    "[W,q,Z];;[W,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;tau2_b1": (100,0,100),
    "[W,q,Z];;[t,g];;pt5;;pt5;;[tracks,allpar,tragam];;tau2_b1": (200,0,200),

    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;tau3_b1": (150,0,150),

    # TAU 1,2,3 B2
    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;tau1_b2": (50,0,25),
    "[t,g];;[W,g,t,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau1_b2": (50,0,50),
    "[t,g];;[W,g,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;tau1_b2": (75,0,75),

    "[g,t];;[W,g,t,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau2_b2": (50,0,50),
    "[W,Z,q];;[W,t,Z,q];;pt1;;pt1;;[tracks,allpar,tragam];;tau2_b2": (40,0,20),
    "[g];;[W,g,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;tau2_b2": (60,0,60),
    "[W,Z,q];;[W,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;tau2_b2": (40,0,40),

    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;tau3_b2": (50,0,50),

    # TAU 3/2, 2/1 *
    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;tau.._b": (100,0,1.1),

    # MULTIPLICITY
    "[g,t,q];;[W,Z,g,t,q];;pt5;;pt5;;[tracks,allpar,tragam];;multiplicity": (50,0,275),

    # C1 B0
    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;c1_b0": (50,0,0.5),

    # C1 B1
    "[W,Z,q];;[W,Z,q];;pt;;pt;;[tracks,allpar,tragam];;c1_b1": (80,0,0.2),
    "[g,t];;[g,t];;pt;;pt;;[tracks,allpar,tragam];;c1_b1": (35,0,0.35),

    # C1 B2
    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;c1_b2": (32,0,0.08),

    # JET PT
    "[W,g,t,Z,q];;[W,g,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;j_pt$": (325,0,7500),

    # JET MASS
    "[W,g,t,Z,q];;[g,t,q];;pt1;;pt1;;[tracks,allpar,tragam];;j_mass_[^m]": (700,0,350),
    "[W];;[Z];;pt1;;pt1;;[tracks,allpar,tragam];;j_mass_[^m]": (350,0,175),
    "[q];;[W,g,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;j_mass_[^m]": (500,0,500),
    "[g,t];;[W,g,t,Z,q];;pt5;;pt5;;[tracks,allpar,tragam];;j_mass_[^m]": (600,0,600),

    # JET MASS MMDT
    "[t];;[W,g,t,Z,q];;pt;;pt;;[tracks,allpar,tragam];;j_mass_mmdt": (100,0,250)
}



### DEFAULT PLOT INFORMATION
#        PLOTNAME          | AXIS TITLE                  | BINS, MIN, MAX |
plotsnames = [
        ('j_pt',           "; jet pT (GeV);",                325, 0, 1500),
        ('j_ptfrac',       "; jet energy fraction;",         60, 0, 2),
        ('j_eta',          "; eta;",                         60, -3, 3),
        ('j_c1_b0',        "; C_{1}^{#beta=0};",             20, 0, 0.5),
        ('j_c1_b1',        "; C_{1}^{#beta=1};",             20, 0, 0.5),
        ('j_c1_b2',        "; C_{1}^{#beta=2};",             20, 0, 0.5),
        ('j_c2_b1',        "; C_{2}^{#beta=1};",             20, 0, 0.5),
        ('j_c2_b2',        "; C_{2}^{#beta=2};",             20, 0, 0.5),
        ('j_d2_b1',        "; D_{2}^{#beta=1};",             20, 0, 0.5),
        ('j_d2_b2',        "; D_{2}^{#beta=2};",             20, 0, 0.5),
        ('j_multiplicity', "; multiplicity;",                40, 0, 200),
        ('j_mass',         "; mass (GeV);",                  80, 0, 300),
        ('j_mass_mmdt',    "; m_{SD}^{#beta=0} (GeV);",      80, 0, 200),
        ('j_mass_sdb2',    "; m_{SD}^{#beta=2} (GeV);",      80, 0, 200),
        ('j_mass_prun',    "; m_{prun} (GeV);",              80, 0, 200),
        ('j_mass_sdm1',    "; m_{SD}^{#beta=1} (GeV);",      80, 0, 200),
        ('j_mass_trim',    "; m_{trim} (GeV);",              80, 0, 200),
        ('j_zlogz',        "; #Sigma z logz;",               28, -6, 1),
        ('j_tau1_b1',      "; N-subjettiness 1, #beta=1;",   100, 0, 150),
        ('j_tau2_b1',      "; N-subjettiness 2, #beta=1;",   100, 0, 150),
        ('j_tau3_b1',      "; N-subjettiness 3, #beta=1;",   100, 0, 150),
        ('j_tau1_b2',      "; N-subjettiness 1, #beta=2;",   100, 0, 150),
        ('j_tau2_b2',      "; N-subjettiness 2, #beta=2;",   100, 0, 150),
        ('j_tau3_b2',      "; N-subjettiness 3, #beta=2;",   100, 0, 150),
        ('j_tau32_b1',     "; N-subjettiness 3/2, #beta=1;", 100, 0, 1),
        ('j_tau21_b1',     "; N-subjettiness 2/1, #beta=1;", 100, 0, 1),
        ('j_tau32_b2',     "; N-subjettiness 3/2, #beta=2;", 100, 0, 1),
        ('j_tau21_b2',     "; N-subjettiness 2/1, #beta=2;", 100, 0, 1)
];

treesToWeight = [
        "t_tracks",
        "t_tragam"
]

branchesToWeight = [
        "j_mass",
        "j_mass_trim",
        "j_mass_sdb2",
        "j_mass_sdm1",
        "j_mass_prun",
        "j_mass_mmdt"
]

def main():
    for n1 in range(0,len(files)) :
        for n2 in range(n1,len(files)) :
            type1,fname1,cat1 = files[n1]
            type2,fname2,cat2 = files[n2]

            if not options.mixCats and cat1 != cat2 : continue

            f1 = ROOT.TFile(fname1);
            f2 = ROOT.TFile(fname2);

            t1_tracks = f1.Get("t_tracks");
            t1_tragam = f1.Get("t_tragam");
            t1_allpar = f1.Get("t_allpar");
            t2_tracks = f2.Get("t_tracks");
            t2_tragam = f2.Get("t_tragam");
            t2_allpar = f2.Get("t_allpar");

            h1_tracks = [];
            h1_tragam = [];
            h1_allpar = [];
            h2_tracks = [];
            h2_tragam = [];
            h2_allpar = [];
            declareHistograms(h1_tracks,"h1_tracks",files[n1],files[n2]);
            declareHistograms(h1_tragam,"h1_tragam",files[n1],files[n2]);
            declareHistograms(h1_allpar,"h1_allpar",files[n1],files[n2]);
            declareHistograms(h2_tracks,"h2_tracks",files[n1],files[n2]);
            declareHistograms(h2_tragam,"h2_tragam",files[n1],files[n2]);
            declareHistograms(h2_allpar,"h2_allpar",files[n1],files[n2]);

            fillHistograms( h1_tracks, t1_tracks,"h1_tracks",("t_tracks" in treesToWeight) );
            fillHistograms( h1_tragam, t1_tragam,"h1_tragam",("t_tragam" in treesToWeight) );
            fillHistograms( h1_allpar, t1_allpar,"h1_allpar",("t_allpar" in treesToWeight) );

            fillHistograms( h2_tracks, t2_tracks,"h2_tracks",("t_tracks" in treesToWeight) );
            fillHistograms( h2_tragam, t2_tragam,"h2_tragam",("t_tragam" in treesToWeight) );
            fillHistograms( h2_allpar, t2_allpar,"h2_allpar",("t_allpar" in treesToWeight) );

            ############################################################

            leg1= ['tracks','tracks and photons','all particles'];
            leg2 = ['%s jet %s'%(type1,cat1),'%s jet %s'%(type2,cat2)]
            nplots = len(plotsnames)
            for i in range(nplots):
                print "Starting i loop... %i"%(i)
                makeCanvas( [h1_tracks[i],h1_tragam[i],h1_allpar[i]],
                        leg1,
                        plotsnames[i][0]+"-%s%s%s-only"%(type1,type1,cat1));
                if n1==n2 : continue
                makeCanvas( [h1_tracks[i],h2_tracks[i]],
                        leg2,
                        plotsnames[i][0]+"-%s%s%sv%s%s%s-tracks"%(type1,type1,cat1,type2,type2,cat2));
                makeCanvas( [h1_tragam[i],h2_tragam[i]],
                        leg2,
                        plotsnames[i][0]+"-%s%s%sv%s%s%s-tragam"%(type1,type1,cat1,type2,type2,cat2));
                makeCanvas( [h1_allpar[i],h2_allpar[i]],
                        leg2,
                        plotsnames[i][0]+"-%s%s%sv%s%s%s-allpar"%(type1,type1,cat1,type2,type2,cat2));

def fillHistograms(hs,t1,tag,applyPtFrac=False):
    #print "Filling histograms %s"%t1.GetName()
    wt = "j_ptfrac[0]" if applyPtFrac else "1";

    for i in xrange(0,len(hs)) :
        #print hs[i].GetName()
        #FilledOnce=False
        for entry in t1.GetListOfBranches():
            #print " - entry = %s"%entry.GetName()
            if "%s_%s"%(entry.GetName(),tag) == hs[i].GetName():
                #if FilledOnce :
                    #%print "CRITICAL ERROR: Filled histogram %s multiple times."%(hs[i].GetName());
                    #%print "Culprit branch name: %s"%(entry.GetName());
                    #if "tau21" in entry.GetName():
                    #    continue
                    #else :
                    #    quit();
                FilledOnce=True
                #print  "--filling!"
                t1.Draw("(%s[0]/%s)>>%s"%(entry.GetName(),(wt if entry.GetName() in branchesToWeight else "1"),hs[i].GetName()));


def declareHistograms(hs,tag,info1,info2):
    print "Declaring histograms %s"%(tag)
    ttype1,tfname1,tcat1 = info1
    ttype2,tfname2,tcat2 = info2

    for varName, axisTitle, nBins, binMin, binMax in plotsnames :
        newNBin=nBins
        newBMin=binMin
        newBMax=binMax

        for regexp in newPlotLimits :
            chkArr=regexp.split(';;');
            passRegexp=ROOT.TString(ttype1).Contains(ROOT.TRegexp(chkArr[0]))
            passRegexp*=ROOT.TString(ttype2).Contains(ROOT.TRegexp(chkArr[1]))
            passRegexp=passRegexp or (ROOT.TString(ttype2).Contains(ROOT.TRegexp(chkArr[0])) and ROOT.TString(ttype1).Contains(chkArr[1]))
            passRegexp*=ROOT.TString(tcat1).Contains(ROOT.TRegexp(chkArr[2]))
            passRegexp*=ROOT.TString(tcat2).Contains(ROOT.TRegexp(chkArr[3]))
            passRegexp*=ROOT.TString(tag.split('_')[1]).Contains(ROOT.TRegexp(chkArr[4]))
            passRegexp*=ROOT.TString(varName).Contains(ROOT.TRegexp(chkArr[5]))

            if passRegexp:
                newNBin,newBMin,newBMax = newPlotLimits[regexp];

        hs.append( ROOT.TH1F("%s_%s"%(varName,tag), axisTitle, newNBin, newBMin, newBMax) );

def makeCanvas(hs,legs,name):
    print "name = ", name;
    colors = [1,2,4,6,7];
    maxval = -999;

    for h in hs:
        if h.Integral() != 0: h.Scale(1/h.Integral());
        h.SetLineWidth(2);
        if h.GetMaximum() > maxval: maxval =  h.GetMaximum()

    leg = ROOT.TLegend(0.2,0.7,0.5,0.9)
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    i = 0;
    for h in hs:
        leg.AddEntry(h,legs[i],"l")
        i+=1;

    c = ROOT.TCanvas("c","c",1000,800);
    hs[0].SetMaximum(maxval*1.5);
    hs[0].Draw("hist");
    i = 0;
    for h in hs:
        h.SetLineColor(colors[i])
        h.Draw("histsames");
        i+=1;
    leg.Draw();

    c.SetLogy(0);
    c.SaveAs(options.outdir+"/"+name+".pdf");
    c.SaveAs(options.outdir+"/"+name+".png");

    if options.makeLogPlots :
        c.SetLogy(1);
        c.SaveAs(options.outdir+"/"+name+"_log.pdf");
        c.SaveAs(options.outdir+"/"+name+"_log.png");

########################################################################################################################
if __name__ == '__main__':
    main();


