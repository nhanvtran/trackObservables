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

parser = OptionParser()

parser.add_option('--inputs',action="store",type="string",dest="inputs",default="./tmp/")
parser.add_option('--output',action="store",type="string",dest="output",default="./")
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
#    c.SaveAs("BackgroundRejection_%s.pdf"%(tree))
#    c.SaveAs("BackgroundRejection_%s.png"%(tree))

# Create plots tracks/allpar ///// tragam/allpar

for pt in [("pt1")] :
    c=ROOT.TCanvas("","",800,800)
    c.cd()

    matHisto=ROOT.TH2F("%smatrix"%(pt), "", procsLen, 0, procsLen, procsLen, 0, procsLen)

    for i in xrange(0,len(procs)) :
        for j in xrange(i+1,len(procs)) :
           allparNum=plotMatrix["allpar"][pt][procs[i]][procs[j]]
           matHisto.Fill(0.5+i, 0.5+j, allparNum/plotMatrix["tracks"][pt][procs[i]][procs[j]])
           matHisto.Fill(0.5+j, 0.5+i, allparNum/plotMatrix["tragam"][pt][procs[i]][procs[j]])

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
    zax.SetRangeUser(0,75)
    zax.SetLabelSize(0.02)


    matHisto.SetTitle("Background Rejection at 50% Signal Efficiency")
    yax.SetTitle("Ratios, %s"%(pt))

    matHisto.Draw("COLZ TEXT")

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

    c.SaveAs("BackgroundRejectionRatios_%s.pdf"%(pt))
    c.SaveAs("BackgroundRejectionRatios_%s.png"%(pt))

for tree,pt,proc1,proc2 in [(t,"pt1",pr1,pr2) for t in trees for pr1 in procs for pr2 in procs]:
    if proc1 == proc2 : continue;
    plotMatrix[tree][pt][proc1][proc2]=plotMatrix["allpar"][pt][proc1][proc2]/plotMatrix[tree][pt][proc1][proc2]
    print "%s \t %s%s \t %s%s \t %.1f" %(tree,proc1,pt,proc2,pt,plotMatrix[tree][pt][proc1][proc2])
