#! /usr/bin/env python

#########################################################
# plotSampleSummary.py                                  #
# Author: E. Coleman, 2016                              #
#                                                       #
# Utility to make a single projection plot with a given #
# input of lines.                                       #
#########################################################

import os
import pprint
import glob
import math
from array import array
import sys
import time
import string
import ROOT

ROOT.gROOT.SetBatch(True)
pp = pprint.PrettyPrinter(indent=4)

import tdrstyle
tdrstyle.setTDRStyle()
#ROOT.gStyle.SetPadTopMargin(0.09);
#ROOT.gStyle.SetPadLeftMargin(0.12);
ROOT.gStyle.SetPadRightMargin(0.05);
#ROOT.gStyle.SetPaintTextFormat("1.1f");

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-i',              action='store',       dest='files',         default="",
        help='file list of inputs ("parton,file,category;...")')
parser.add_option('-n',              action='store',       dest='outname',        default="",
        help='where to store plotting output')
parser.add_option('-o',              action='store',       dest='outdir',        default="./plots/",
        help='where to store plotting output')
parser.add_option('--logPlots',      action='store_false', dest='makeLogPlots',  default=True,
        help='make log plots')
parser.add_option('--basedir',       action='store',       dest='base',          default='../testSamples/',
        help='location of input root files to plot')
parser.add_option('--ana',           action='store',       dest='anasub',
        default="r0_h0_e0,r05_h05_e005,r05_h01_e005,r05_h01_e005_t,r05_h002_e005_t",
        help='which anaSubstructures to use (csv)')
parser.add_option('--pts',           action='store',       dest='ptList',        default="pt1,pt5",
        help='which pts to use (csv)')
parser.add_option('--sigs',          action='store',       dest='sigList',       default="W,Z,t,q,g",
        help='which signals to use (csv)')
parser.add_option('--trees',         action='store',       dest='treeList',      default="allpar,tragam,tracks",
        help='which trees to use (csv)')
parser.add_option('--vars',          action='store',       dest='vars',          default="j_mass_mmdt",
        help='which variables to plot')
parser.add_option('--lines',         action='store',       dest='lines',         default="",
        help='line command to plot')

(options, args) = parser.parse_args()

if options.lines=="" :
    print "ERROR: No lines specified"
    exit(0)

# useful tags (HARDCODED)
treeNames = {
            'tracks': 'Tracks',
            'tragam': 'Tracks+#gamma',
            'allpar': 'All par.'
            }
anaNames={"r0_h0_e0"         : "Perfect" ,
        "r05_h05_e005"       : "HCAL0.05",
        "r05_h01_e005"       : "HCAL0.01",
        "r05_h01_e005_t"     : "F1",
        "r05_h002_e001_t"    : "F2",
        "r05_h01_e005_t220"  : "F1",
        "r05_h002_e001_t220" : "F2",
        "r05_h01_e005_t500"  : "F1",
        "r05_h002_e001_t500" : "F2",
        "r05_h002_e005_t500" : "Extreme-gran.",
        "r05_h005_e005_t500" : "High-gran.",
         "r1_h022_e0175_t220": "CMS-like",
         "r1_h022_e0175_t220_q10": "CMS-like (+PU)",
         "r1_h022_e050_t110" : "CMS-like (EB)",
         "r1_h022_e0175_t500": "CMS-like (500 deg.)",
         "r1_h022_e005_t110" : "CMS-like (110 deg, E0.005)",
         "r1_h022_e0175_t110": "CMS-like (110 deg.)"}
ptNames = {
            'pt1': 'p_{T} 1 TeV',
            'pt5': 'p_{T} 5 TeV'
            }
ptTags = {
            'pt1': 'processed-pythia82-lhc13',
            'pt5': 'processed-pythia82-fcc100'
            }

# helpers
vars   =options.vars.split(',')
anaSubs=options.anasub.split(',')
pts    =options.ptList.split(',')
sigs   =options.sigList.split(',')
trees  =options.treeList.split(',')
base   =options.base
files=[]

#####################################################
# function to help with automation: replace * in input
# lines with all possible options
def replaceWildCards(map,searchChar="*",copyLines=False):
    tPlotMap=[]
    plotAgain=False
    # scan through lines
    for plot in map :
        tPlot = [[]]
        for line in plot :
            # if no wildcard, just add line
            i = line.index(searchChar)
            replMap=[]
            if   i==0: replMap=trees
            elif i==1: replMap=pts
            elif i==2: replMap=sigs
            elif i==3: replMap=anaSubs

            if searchChar not in line or len(line)==1 :
                for cPlot in plotMap :
                    cPlot += line
                continue
            # else replace the * with the relevant array

            # perform the replacement
            irepl=0
            for repl in replMap :
                tdata=line[:]
                tdata[i]=repl
                if copyLines :
                    for inew in range(len(tPlot),irepl+1) :
                        tPlot += [[]]
                    tPlot[irepl] += [tdata]
                else : tPlot[0]+=[tdata]
                irepl+=1
        tPlotMap += tPlot

    if searchChar in [datum for thing in tPlot for line in thing for datum in line]:
        plotAgain=True

    # recurse if still wildcards in the array
    if plotAgain: tPlotMap=replaceWildCards(tPlotMap,searchChar,copyLines)
    return tPlotMap
#####################################################

# prep files
for ana,sig,pt in [(a,b,c)
        for a in anaSubs
        for b in sigs
        for c in pts] :
    if options.files!="" : break
    anasub="-"+ana if options.anasub != "" else ""
    files+=[(sig,'%s/%s-%s%s-%s-50k%s.root'%(base,ptTags[pt],sig,sig,pt,anasub),pt)]

if options.files!="" :
    files=[(x.split(',')[0],x.split(',')[1],x.split(',')[2]) for x in options.files.split(';')];


def getFile(pt,sig,ana) :
    for tsig,name,tpt in files :
        if sig != sig : continue
        if pt  != tpt : continue
        if sig in name and pt in name and ana in name :
            return ROOT.TFile(name,"READ")
    return None


# var;;[tree1,pt1,sig1,ana1];;[tree2,pt2,sig2,ana2];;etc.
# can use * to reproduce for all
plotMap=[[x.split(',') for x in options.lines.split(';')]]
plotMap=replaceWildCards(plotMap,searchChar="**",copyLines=True)
plotMap=replaceWildCards(plotMap)

pp.pprint(plotMap)


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
        ('j_pt',           "; jet pT (GeV);",                325, 0, 7500),
        ('j_ptfrac',       "; jet energy fraction;",         60, 0, 2),
        ('j_eta',          "; eta;",                         60, -3, 3),
        ('j_c1_b0',        "; C_{1}^{#beta=0};",             20, 0, 0.5),
        ('j_c1_b1',        "; C_{1}^{#beta=1};",             20, 0, 0.5),
        ('j_c1_b2',        "; C_{1}^{#beta=2};",             20, 0, 0.5),
        ('j_c2_b1',        "; C_{2}^{#beta=1};",             20, 0, 0.5),
        ('j_c2_b2',        "; C_{2}^{#beta=2};",             20, 0, 0.5),
        ('j_d2_b1',        "; D_{2}^{#beta=1};",             70, 0, 8.0),
        ('j_d2_b2',        "; D_{2}^{#beta=2};",             200, 0, 50.0),
        ('j_multiplicity', "; multiplicity;",                40, 0, 200),
        ('j_mass',         "; mass (GeV);",                  80, 0, 300),
        ('j_mass_mmdt',    "; m_{SD}^{#beta=0} (GeV);",      80, 0, 125),
        ('j_mass_sdb2',    "; m_{SD}^{#beta=2} (GeV);",      80, 0, 200),
        ('j_mass_prun',    "; m_{prun} (GeV);",              80, 0, 200),
        ('j_mass_sdm1',    "; m_{SD}^{#beta=1} (GeV);",      80, 0, 200),
        ('j_mass_trim',    "; m_{trim} (GeV);",              80, 0, 200),
        ('j_zlogz',        "; #Sigma z logz;",               28, -6, 1),
        ('j_tau1_b1',      "; #tau_{1}^{#beta=1};",   100, 0, 150),
        ('j_tau2_b1',      "; #tau_{2}^{#beta=1};",   100, 0, 150),
        ('j_tau3_b1',      "; #tau_{3}^{#beta=1};",   100, 0, 150),
        ('j_tau1_b2',      "; #tau_{1}^{#beta=2};",   100, 0, 150),
        ('j_tau2_b2',      "; #tau_{2}^{#beta=2};",   100, 0, 150),
        ('j_tau3_b2',      "; #tau_{3}^{#beta=2};",   100, 0, 150),
        ('j_tau32_b1',     "; #tau_{32}^{#beta=1};", 80, 0, 1),
        ('j_tau21_b1',     "; #tau_{21}^{#beta=1};", 80, 0, 1),
        ('j_tau32_b2',     "; #tau_{32}^{#beta=2};", 80, 0, 1),
        ('j_tau21_b2',     "; #tau_{21}^{#beta=2};", 80, 0, 1)
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
    arrOfNameables= [ i
            for i, x in enumerate(options.lines.split(';')[0].split(',')) if x=="**"]

    for plotLines,var in [(a,b)
            for a in plotMap
            for b in vars]:
        # [tree1,pt1,sig1,ana1];[tree2,pt2,sig2,ana2];etc.

        tname=""
        hists=[]
        leg=[]

        for i in range(len(arrOfNameables)) :
            if i in arrOfNameables : tname += "_"+plotLines[0][i]

        sameMap=[""]*4
        for tree,pt,sig,ana in [tuple(x) for x in plotLines]:
            if sameMap[0] == "" :
                sameMap=[sig,treeNames[tree],ptNames[pt],anaNames[ana]] # change order because of legend format
                continue
            if sameMap[0] != sig and sameMap[0] is not None :
                sameMap[0] = None
                continue
            if sameMap[1] != treeNames[tree] and sameMap[1] is not None :
                sameMap[1] = None
                continue
            if sameMap[2] != ptNames[pt] and sameMap[2] is not None:
                sameMap[2] = None
                continue
            if sameMap[3] != anaNames[ana] and sameMap[3] is not None:
                sameMap[3] = None
                continue


        for tree,pt,sig,ana in [tuple(x) for x in plotLines]:

            fIn=getFile(pt,sig,ana)
            ttree=fIn.Get("t_"+tree)
            thist=declareHistogram(pt,sig,tree,var)

            legTitle=[]
            if sig not in sameMap :  legTitle+=[sig]
            if treeNames[tree] not in sameMap : legTitle+=[treeNames[tree]]
            if ptNames[pt] not in sameMap :   legTitle+=[ptNames[pt]]
            if anaNames[ana] not in sameMap :  legTitle+=[anaNames[ana]]
            leg+=[', '.join(legTitle)]

            wt="j_ptfrac[0]" if ("t_"+tree in treesToWeight and var in branchesToWeight) else "1"
            ttree.Draw("(%s/%s)>>%s"%(var,wt,thist.GetName()))
            hists+=[thist.Clone(thist.GetName())]
            for hist in hists : hist.SetDirectory(0)
            print ""

        pp.pprint(hists)

        tname='_'.join([x for x in sameMap if x is not None])
        tname=tname.replace('#','')
        tname=tname.replace(' ','')
        tname=tname.replace('+','')
        tname=tname.replace('{','')
        tname=tname.replace('}','')
        tname=tname.replace('p_T','pt')
        tname=tname.replace('.','p')
        makeCanvas(hists,leg,"SummaryPlot_%s_%s_%s"%(options.outname,var,tname),otherTextMap=sameMap)
        del hists
        del leg
        ROOT.gROOT.CloseFiles()

def declareHistogram(pt,tree,sig,var):
    for varName, axisTitle, nBins, binMin, binMax in plotsnames :
        if varName != var : continue
        newNBin=nBins
        newBMin=binMin
        newBMax=binMax

        #for regexp in newPlotLimits :
        #    chkArr=regexp.split(';;');
        #    passRegexp =ROOT.TString(sig    ).Contains(ROOT.TRegexp(chkArr[0]))
        #    passRegexp*=ROOT.TString(pt     ).Contains(ROOT.TRegexp(chkArr[2]))
        #    passRegexp*=ROOT.TString(tree   ).Contains(ROOT.TRegexp(chkArr[4]))
        #    passRegexp*=ROOT.TString(varName).Contains(ROOT.TRegexp(chkArr[5]))

        #    if passRegexp:
        #        newNBin,newBMin,newBMax = newPlotLimits[regexp];

        hist=ROOT.TH1F("%s_%s_%s_%s"%(varName,sig,pt,tree), axisTitle, newNBin, newBMin, newBMax)
        return hist
    return None

def makeCanvas(hs,legs,name,otherTextMap=[]):
    print "name = ", name;
    colors = [ROOT.kRed+2,ROOT.kBlue+2,ROOT.kRed,ROOT.kBlue+1,ROOT.kOrange+7,ROOT.kAzure+7] # [1,2,4,6,7,8,9,10,11];
    #colors = [ROOT.kRed,ROOT.kOrange+7,ROOT.kBlue+1,ROOT.kAzure+7,ROOT.kGreen+2,ROOT.kGreen] # [1,2,4,6,7,8,9,10,11];
    maxval = -999;

    for h in hs:
        if h.Integral() != 0: h.Scale(1/h.Integral());
        h.SetLineWidth(2);
        if h.GetMaximum() > maxval: maxval =  h.GetMaximum()

    leg = ROOT.TLegend(0.2,0.5,0.5,0.7)
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    i = 0;
    for h in hs:
        leg.AddEntry(h,legs[i],"l")
        i+=1;

    c = ROOT.TCanvas("c","c",800,800);
    hs[0].SetMaximum(maxval*1.5);
    hs[0].Draw("c hist");
    i = 0;
    for h in hs:
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        h.Draw("c histsames");
        i+=1;
    leg.Draw();

    otherText=""
    for text in otherTextMap :
        if text is None : continue
        if otherText=="" : otherText=text
        else : otherText="#splitline{%s}{%s}"%(otherText,text + (" samples" if len(text)==1 else ""))
    otherLatex=ROOT.TLatex()
    otherLatex.SetTextSize(0.04)
    otherLatex.DrawLatexNDC(0.225,0.927 - 0.06*len([x for x in otherTextMap if x is not None]),otherText)

    c.SetLogy(0);

    c.SaveAs(options.outdir+"/"+name+".pdf");
    c.SaveAs(options.outdir+"/"+name+".png");
    c.SaveAs(options.outdir+"/"+name+".C");

    if options.makeLogPlots :
        c.SetLogy(1);
        c.SaveAs(options.outdir+"/"+name+"_log.pdf");
        c.SaveAs(options.outdir+"/"+name+"_log.png");
        c.SaveAs(options.outdir+"/"+name+"_log.C");
    del c
    del leg

########################################################################################################################
if __name__ == '__main__':
    main();


