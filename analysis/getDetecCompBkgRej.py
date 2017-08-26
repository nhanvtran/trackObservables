#! /usr/bin/env python
import pprint
import os
import glob
import math
from array import array
import sys
import time

############################################
#            Job steering                  #
############################################
from optparse import OptionParser
import ROOT

ROOT.gROOT.SetBatch(True)

parser = OptionParser()

parser.add_option('--inputs',action="store",type="string",dest="inputs",default="./tmp/")
parser.add_option('--output',action="store",type="string",dest="output",default="./")
parser.add_option('--training',action="store",type="string",dest="training",default="allcut")
parser.add_option('--tree',action="store",type="string",dest="tree",default="allpar")
parser.add_option('-n', '--name',action="store",type="string",dest="outname",default="DetecComp_bkgRejPwr_at_Sig0p5")

(options, args) = parser.parse_args()

onlyfiles=glob.glob(options.inputs)
fout = open('%s/%s.txt'%(options.output,options.outname),'w')

treeNames={ "allpar" : "All particles",
            "tragam" : "Tracks+#gamma",
            "tracks" : "Tracks only" }
treename = treeNames[options.tree]

pts=['pt1','pt5']
procs=["W","Z","q","g","t"]
procsLen = len(procs)

plotMatrix={}
ptNames={ "pt1" : "p_{T} 1 TeV",
        "pt5"   : "p_{T} 5 TeV"}
anas = ["CMS-like", "F1", "F2"]
anaNames={"r0_h0_e0"       :   "Perfect" ,
        "r05_h01_e005_t220":   "F1",
        "r05_h002_e001_t220":  "F2",
         "r1_h022_e0175_t220": "CMS-like"}
trainingNames={ "all" : "All observables",
        "shapesonly"  : "Shape observables",
        "massonly"    : "Mass observables",
        "allcut"         : "All observables",
        "shapesonlycut"  : "Shapes observables",
        "massonlycut"    : "Mass observables" }
trainingName="" if options.training not in trainingNames else trainingNames[options.training]

# Prepare arrays
for ana in anas:
    plotMatrix[ana]={}
    for pt in pts:
        plotMatrix[ana][pt]={}
        for proc in procs:
            plotMatrix[ana][pt][proc]={}

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(onlyfiles)
pp.pprint(plotMatrix)

# Fill arrays and make output txt file
for tfname in onlyfiles :
    print tfname
    if options.tree not in tfname : continue

    anaName = tfname.split('/')[-4]
    anaName = anaNames[anaName]

    fname = tfname.split('/')[-1]
    fIn=ROOT.TFile(tfname)

    MVAh=fIn.Get("Method_BDT/BDTG/MVA_BDTG_effBvsS")
    halfBin=MVAh.GetXaxis().FindBin(0.60)
    cat1=fname.split("_")[2]
    cat2=fname.split("_")[3]
    treename=fname.split("_")[-1].replace('.root','')
    print "\t - %s %s %s %s %s"%(anaName,treename,cat1.split('-')[1],cat1[0],cat2[0])

    if MVAh.GetBinContent(halfBin) == 0 :
        print "ERROR: bin content is 0"
        continue
    bkgRej=1/MVAh.GetBinContent(halfBin)
    plotMatrix[anaName][cat1.split('-')[1]][cat1[0]][cat2[0]]=bkgRej
    fout.write("%s \t %s \t %s \t %.3f\n"%(anaName, cat1, cat2, bkgRej))

    fIn.Close()

fout.close()

exit

pp.pprint(plotMatrix)

# Create plots F1/CMS ///// F2/CMS

for pt in ["pt1","pt5"] :
    c=ROOT.TCanvas("","",800,800)
    c.cd()

    matHisto=ROOT.TH2F("%smatrix"%(pt), "", procsLen, 0, procsLen, procsLen, 0, procsLen)

    for i in xrange(0,len(procs)) :
        for j in xrange(i+1,len(procs)) :
           print procs[i],procs[j]
           try :
               CMSnum=plotMatrix["CMS-like"][pt][procs[i]][procs[j]]
               matHisto.Fill(0.5+i, 0.5+j, plotMatrix["F1"][pt][procs[i]][procs[j]]/CMSnum)
               matHisto.Fill(0.5+j, 0.5+i, plotMatrix["F2"][pt][procs[i]][procs[j]]/CMSnum)
           except KeyError :
               print "WARNING: Need to switch procs!"
               CMSnum=plotMatrix["CMS-like"][pt][procs[j]][procs[i]]
               matHisto.Fill(0.5+i, 0.5+j, plotMatrix["F1"][pt][procs[j]][procs[i]]/CMSnum)
               matHisto.Fill(0.5+j, 0.5+i, plotMatrix["F2"][pt][procs[j]][procs[i]]/CMSnum)


    xax=matHisto.GetXaxis()
    xax.SetTitle("")
    yax=matHisto.GetYaxis()
    yax.SetTitle("")
    for i in xrange(0,procsLen) :
        proc = procs[i]
        bin_index_x = xax.FindBin(0.5+i)
        bin_index_y = yax.FindBin(0.5+i)
        label = "%s%s"%(proc,proc)
        xax.SetBinLabel(bin_index_x,label)
        yax.SetBinLabel(bin_index_y,label)

    ROOT.gStyle.SetPaintTextFormat(".1f")
    ROOT.gStyle.SetOptStat(0)

    zax=matHisto.GetZaxis()
    zax.SetRangeUser(0.99,75)
    zax.SetLabelSize(0.02)


    matHisto.SetTitle("Background Rejection at 60% Signal Efficiency")
    yax.SetTitle("%s Trainings, %s, %s"%(trainingName,ptNames[pt],treeNames[options.tree]))
    yax.SetTitleOffset(1.2)

    matHisto.Draw("COLZ TEXT")
    c.SetLogz()

    sepline=ROOT.TLine()
    sepline.SetLineColor(ROOT.kRed)
    sepline.SetLineWidth(2)
    sepline.DrawLineNDC(0.1,0.1,0.9,0.9)

    pt1text=ROOT.TLatex()
    pt1text.SetTextColor(ROOT.kRed)
    pt1text.SetTextAngle(45)
    pt1text.SetTextSize(0.03)
    pt1text.DrawLatexNDC(.435,.446,"F1 / CMS-like")

    pt5text=ROOT.TLatex()
    pt5text.SetTextColor(ROOT.kRed)
    pt5text.SetTextAngle(45)
    pt5text.SetTextSize(0.03)
    pt5text.DrawLatexNDC(.458,.42,"F2 / CMS-like")

    c.SaveAs("%s/DetecComp_BackgroundRejectionRatios60_%s_%s_%s.pdf"%(options.output,pt,options.tree,options.training))
    c.SaveAs("%s/DetecComp_BackgroundRejectionRatios60_%s_%s_%s.png"%(options.output,pt,options.tree,options.training))
    c.SaveAs("%s/DetecComp_BackgroundRejectionRatios60_%s_%s_%s.C"  %(options.output,pt,options.tree,options.training))
