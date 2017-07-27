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
parser.add_option('--ana'   ,action="store",type="string",dest="ana",default="")
parser.add_option('--training',action="store",type="string",dest="training",default="")
parser.add_option('-n', '--name',action="store",type="string",dest="outname",default="bkgRejPwr_at_Sig0p5")

(options, args) = parser.parse_args()

onlyfiles=os.listdir(options.inputs)
fout = open('%s/%s.txt'%(options.output,options.outname),'w')

trees=["tragam","tracks","allpar"]
pts=['pt1','pt5']
procs=["W","Z","q","t","g"]
procsLen = len(procs)

plotMatrix={}
maxVal=0

ptNames={ "pt1" : "p_{T} 1 TeV",
        "pt5"   : "p_{T} 5 TeV"}
anaNames={"r0_h0_e0"       :   "Perfect" ,
        "r05_h05_e005"     :   "HCAL0.05",
        "r05_h01_e005"     :   "HCAL0.01",
        "r05_h01_e005_t":   "F1",
        "r05_h002_e001_t":  "F2",
        "r05_h01_e005_t500":   "F1",
        "r05_h002_e001_t500":  "F2",
        "r05_h002_e005_t500":  "Extreme-gran.",
        "r05_h005_e005_t500":  "High-gran.",
         "r1_h022_e050_t110":  "CMS-like (EB)",
         "r1_h022_e0175_t110": "CMS-like"}
trainingNames={ "all" : "All Variables",
        "shapesonly"  : "Shape Variables",
        "massonly"    : "Mass Variables"}
trainingName="" if options.training not in trainingNames else trainingNames[options.training]
anaName = None if options.ana not in anaNames else anaNames[options.ana]

# Prepare arrays
for tree in trees:
    plotMatrix[tree]={}
    for pt in pts:
        plotMatrix[tree][pt]={}
        for proc in procs:
            plotMatrix[tree][pt][proc]={}

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(onlyfiles)
pp.pprint(plotMatrix)

# Fill arrays and make output txt file
for fname in onlyfiles :
    print fname
    fIn=ROOT.TFile("%s/%s"%(options.inputs,fname))

    MVAh=fIn.Get("Method_BDT/BDTG/MVA_BDTG_effBvsS")
    halfBin=MVAh.GetXaxis().FindBin(0.50)

    cat1=fname.split("_")[2]
    cat2=fname.split("_")[3]
    treename=fname.split("_")[-1].replace('.root','')
    print "\t - %s %s %s %s"%(treename,cat1.split('-')[1],cat1[0],cat2[0])

    if MVAh.GetBinContent(halfBin) == 0 :
        print "ERROR: bin content is 0"
        continue
    bkgRej=1/MVAh.GetBinContent(halfBin)
    plotMatrix[treename][cat1.split('-')[1]][cat1[0]][cat2[0]]=bkgRej
    #plotMatrix[treename][cat1.split('-')[1]][cat2[0]][cat1[0]]=bkgRej
    fout.write("%s \t %s \t %s \t %.3f\n"%(cat1, cat2, treename, bkgRej))

    fIn.Close()

fout.close()

pp.pprint(plotMatrix)

# Create plots pt1 / pt5

#for tree in [(tree) for tree in trees] :
#    c=ROOT.TCanvas("","",800,800)
#    c.cd()
#
#    matHisto=ROOT.TH2F("%smatrix"%(tree), "", procsLen, 0, procsLen, procsLen, 0, procsLen)
#
#    for i in xrange(0,len(procs)) :
#        for j in xrange(i+1,len(procs)) :
#           matHisto.Fill(0.5+i, 0.5+j, plotMatrix[tree]["pt1"][procs[i]][procs[j]])
#           matHisto.Fill(0.5+j, 0.5+i, plotMatrix[tree]["pt5"][procs[i]][procs[j]])
#
#    xax=matHisto.GetXaxis()
#    xax.SetTitle("")
#    yax=matHisto.GetYaxis()
#    yax.SetTitle("")
#    for i in xrange(0,procsLen) :
#        proc = procs[i]
#        bin_index_x = xax.FindBin(0.5+i)
#        bin_index_y = yax.FindBin(0.5+i)
#        label = "%s%s"%(proc,proc)
#        xax.SetBinLabel(bin_index_x,label)
#        yax.SetBinLabel(bin_index_y,label)
#
#    ROOT.gStyle.SetPaintTextFormat(".1f")
#    ROOT.gStyle.SetOptStat(0)
#
#    zax=matHisto.GetZaxis()
#    zax.SetRangeUser(0,3200)
#    zax.SetLabelSize(0.02)
#
#
#    matHisto.SetTitle("Background Rejection at 50% Signal Efficiency")
#    yax.SetTitle("%s"%(tree))
#
#    matHisto.Draw("COLZ TEXT")
#
#    sepline=ROOT.TLine()
#    sepline.SetLineColor(ROOT.kRed)
#    sepline.SetLineWidth(2)
#    sepline.DrawLineNDC(0.1,0.1,0.9,0.9)
#
#    pt1text=ROOT.TLatex()
#    pt1text.SetTextColor(ROOT.kRed)
#    pt1text.SetTextAngle(45)
#    pt1text.SetTextSize(0.025)
#    pt1text.DrawLatexNDC(.44,.48,"p_{T} 1 TeV")
#
#    pt5text=ROOT.TLatex()
#    pt5text.SetTextColor(ROOT.kRed)
#    pt5text.SetTextAngle(45)
#    pt5text.SetTextSize(0.025)
#    pt5text.DrawLatexNDC(.51,.45,"p_{T} 5 TeV")
#
#    c.SaveAs("%s/BackgroundRejection_%s.pdf"%(options.output,tree))
#    c.SaveAs("%s/BackgroundRejection_%s.png"%(options.output,tree))

# Create plots tracks/allpar ///// tragam/allpar

for pt in ["pt1","pt5"] :
    c=ROOT.TCanvas("","",800,800)
    c.cd()

    matHisto=ROOT.TH2F("%smatrix"%(pt), "", procsLen, 0, procsLen, procsLen, 0, procsLen)

    for i in xrange(0,len(procs)) :
        for j in xrange(i+1,len(procs)) :
           print procs[i],procs[j]
           try :
               allparNum=plotMatrix["allpar"][pt][procs[i]][procs[j]]
               matHisto.Fill(0.5+i, 0.5+j, allparNum/plotMatrix["tracks"][pt][procs[i]][procs[j]])
               matHisto.Fill(0.5+j, 0.5+i, allparNum/plotMatrix["tragam"][pt][procs[i]][procs[j]])
           except KeyError :
               print "WARNING: Need to switch procs!"
               allparNum=plotMatrix["allpar"][pt][procs[j]][procs[i]]
               matHisto.Fill(0.5+i, 0.5+j, allparNum/plotMatrix["tracks"][pt][procs[j]][procs[i]])
               matHisto.Fill(0.5+j, 0.5+i, allparNum/plotMatrix["tragam"][pt][procs[j]][procs[i]])


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
    zax.SetRangeUser(1,75)
    zax.SetLabelSize(0.02)


    matHisto.SetTitle("Background Rejection at 50% Signal Efficiency")
    yax.SetTitle("%s Trainings, %s%s"%(trainingName,ptNames[pt],("" if anaName is None else ", "+anaName)))
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
    pt1text.SetTextSize(0.025)
    pt1text.DrawLatexNDC(.435,.455,"all par. / tracks")

    pt5text=ROOT.TLatex()
    pt5text.SetTextColor(ROOT.kRed)
    pt5text.SetTextAngle(45)
    pt5text.SetTextSize(0.025)
    pt5text.DrawLatexNDC(.45,.42,"all par. / tracks+#gamma")

    c.SaveAs("%s/BackgroundRejectionRatios_%s.pdf"%(options.output,pt))
    c.SaveAs("%s/BackgroundRejectionRatios_%s.png"%(options.output,pt))

for tree,pt,proc1,proc2 in [(t,pt,pr1,pr2)
        for t in trees
        for pt in pts
        for pr1 in procs
        for pr2 in procs]:
    if proc1 == proc2 : continue;
    print "%s \t %s%s \t %s%s \t %.1f" %(tree,proc1,pt,proc2,pt,plotMatrix[tree][pt][proc1][proc2])
    plotMatrix[tree][pt][proc1][proc2]=plotMatrix["allpar"][pt][proc1][proc2]/plotMatrix[tree][pt][proc1][proc2]
