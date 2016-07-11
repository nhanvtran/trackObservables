#! /usr/bin/env python
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

for fname in onlyfiles :
    print fname
    fIn=ROOT.TFile("%s/%s"%(options.inputs,fname))

    MVAh=fIn.Get("Method_BDT/BDTG/MVA_BDTG_effBvsS")
    halfBin=MVAh.GetXaxis().FindBin(0.50)

    cat1=fname.split("_")[2]
    cat2=fname.split("_")[3]
    treename=fname.split("_")[-1].replace('.root','')

    fout.write("%s \t %s \t %s \t %.3f\n"%(cat1, cat2, treename, 1/MVAh.GetBinContent(halfBin)))

    fIn.Close()

fout.close()
