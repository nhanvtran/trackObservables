#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import pprint
import ROOT

ROOT.gROOT.SetBatch(False)
pp = pprint.PrettyPrinter(indent=2)


def GetEnvelope(aloTH1) :
    nPoints = aloTH1[0].GetXaxis().GetNbins()+1
    x=ROOT.TVector(nPoints)
    y=ROOT.TVector(nPoints)
    xerr=ROOT.TVector(nPoints)
    yerr=ROOT.TVector(nPoints)

    for iBin in range(0,nPoints) :
        maxVal=0
        minVal=1
        for histo in aloTH1 :
            maxVal = max(maxVal,histo.GetBinContent(iBin))
            minVal = min(minVal,histo.GetBinContent(iBin))

        y[iBin]    = (maxVal+minVal)/2
        yerr[iBin] = (maxVal-minVal)/2
        x[iBin]    = aloTH1[0].GetXaxis().GetBinCenter(iBin)
        xerr[iBin] = 0

    gr = ROOT.TGraphErrors(x,y,xerr,yerr)
    return gr


############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--indir',action="store",type="string",dest="indir",default="./output/")
parser.add_option('--outdir',action="store",type="string",dest="outdir",default="./output/")
parser.add_option('--pts',action="store",type="string",dest="pts",default="pt1")
parser.add_option('--trainings',action="store",type="string",dest="trainings",default="all")
parser.add_option('--anasubs',action="store",type="string",dest="anasubs",default="r1_h022_e0175_t220,r05_h01_e005_t220,r05_h002_e001_t220")
parser.add_option('--infos',action="store",type="string",dest="infos",default="allpar,tragam,tracks")
parser.add_option('--procs',action="store",type="string",dest="procs",default="W,q;W,Z;q,g")

(options, args) = parser.parse_args()

anaSubList   = options.anasubs.split(',')
trainingList = options.trainings.split(',')
ptList       = options.pts.split(',')
procList     = [x.split(',') for x in options.procs.split(';')]
infoList     = options.infos.split(',')

infoNames = {
        "tracks": "Tracks",
        "tragam": "Tracks+#gamma",
        "allpar": "All par."
        }
trainingNames = {
        "all" : "All vars."
        }
ptNames = {
        "pt5" : "p_{T} 5 TeV",
        "pt1" : "p_{T} 1 TeV"
        }

canvas = ROOT.TCanvas("cROC","cROC",700,700)
canvas.SetLogy()

ROOT.gStyle.SetOptTitle(False)
ROOT.gStyle.SetPalette(ROOT.kSolar)

colorList = [ROOT.kRed+3,ROOT.kRed+1,ROOT.kOrange+7]

for procs in procList :

    sg=procs[0]
    bg=procs[1]

    iColor=0
    canvas.Clear()

    mg = ROOT.TMultiGraph()

    for info in infoList:
        cHistos=[]

        for pt,anasub,training in [(a,b,c)
                for a in ptList
                for b in anaSubList
                for c in trainingList] :

            fLoc = options.indir+"/"+anasub+"/eosrootfiles/"+training
            fLoc += "/MVA_bdtg_"+sg+sg+'-'+pt+'_'
            fLoc += bg+bg+'-'+pt+'_t_'+info+'.root'
            fIn = ROOT.TFile(fLoc)

            fIn.Print()
            ROC = fIn.Get("Method_BDT/BDTG/MVA_BDTG_effBvsS").Clone()
            ROC.SetDirectory(0)
            cHistos+=[ROC]

            fIn.Close()

        print cHistos

        curGraph=GetEnvelope(cHistos).Clone()
        curGraph.SetFillColor(colorList[iColor])
        curGraph.SetTitle(infoNames[info])
        #curGraph.Draw("a4"+(" SAME" if not iColor == 0 else ""))
        iColor += 1
        mg.Add(curGraph,"4")

    mg.Draw("A")
    mg.GetXaxis().SetRangeUser(0,1);
    mg.GetYaxis().SetRangeUser(5*10e-06,1);
    mg.GetXaxis().SetTitle("Signal Efficiency")
    mg.GetYaxis().SetTitle("Background Efficiency")
    mg.GetXaxis().SetTitleSize(0.04);
    mg.GetYaxis().SetTitleSize(0.04);
    mg.GetXaxis().SetTitleOffset(1.05);
    mg.GetYaxis().SetTitleOffset(1.2);
    mg.GetXaxis().SetLabelSize(0.03);
    mg.GetYaxis().SetLabelSize(0.03);
    canvas.SetTickx(1);
    canvas.SetTicky(1);
    canvas.SetRightMargin(0.05);
    canvas.SetTopMargin(0.05);
    leg=canvas.BuildLegend(0.15,0.68,0.47,0.88,"")
    leg.SetBorderSize(0);
    leg.SetTextSize(0.035);
    leg.SetTextFont(42);
    leg.SetLineColor(1);
    leg.SetLineStyle(1);
    leg.SetLineWidth(1);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);

    extra = "#splitline{"+trainingNames[training]+"}{"+ptNames[pt]+"}";
    extra = "#splitline{"+sg+" vs. "+bg+"}{"+extra+"}";

    latex = ROOT.TLatex();
    latex.SetTextSize(0.04);
    latex.DrawLatexNDC(0.725,0.25,extra);

    canvas.Modified()
    canvas.Update()

    raw_input()
    canvas.SaveAs("ROCEnvelope_"+sg+bg+"_"+pt+".png")
    canvas.SaveAs("ROCEnvelope_"+sg+bg+"_"+pt+".pdf")
