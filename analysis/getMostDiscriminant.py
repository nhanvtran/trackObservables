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

parser = OptionParser()

parser.add_option('--inputs',    action="store",type="string",dest="inputs",default="./output/")
parser.add_option('--anasubs',   action="store",type="string",dest="anasubs",default="r05_h02_e005,r05_h01_e005,r05_h01_e005_t")
parser.add_option('--trainings', action="store",type="string",dest="trainings",default="all,shapesonly,massonly")
parser.add_option('--sigs',      action="store",type="string",dest="sigs", default="tt,WW,ZZ,gg,qq")
parser.add_option('--bkgs',      action="store",type="string",dest="bkgs", default="tt,WW,ZZ,gg,qq")
parser.add_option('--trees',     action="store",type="string",dest="trees",default="tracks,tragam,allpar")
parser.add_option('--output',    action="store",type="string",dest="output",default="./")
parser.add_option('-n',          action="store",type="string",dest="outname",default="tot")

(options, args) = parser.parse_args()

varCounts={}
varCounts['pt1']={}
varCounts['pt5']={}

if __name__ == '__main__':

    anaSubList=options.anasubs.split(',')
    trainList =options.trainings.split(',')
    treeList  =options.trees.split(',')
    sigList   =options.sigs.split(',')
    bkgList   =options.bkgs.split(',')

    for ana,train,sig,bkg,tree in [(a,b,c,d,e)
            for a in anaSubList
            for b in trainList
            for c in sigList
            for d in bkgList
            for e in treeList]:
        if sig == bkg : continue
        foundTable = False;
        numBars = 0;
        with open("%s/%s/%s/%s-%s-t_%s.txt"%(options.inputs,ana,train,sig,bkg,tree), 'r') as fin :
            empty=False
            for i, line in enumerate(fin):
                if i == 0 and line == '' : empty=True
                else : break
            if empty : continue


            rank=0
            curPt="pt1"
            for line in fin:
                if line == "\n" or line=="" : continue
                if "pt 1" in line:
                    curPt="pt1"
                    continue
                if "pt 5" in line:
                    curPt="pt5"
                    continue

                infoLine=line
                if 'Rank' in infoLine :
                    foundTable=True
                    continue
                if foundTable :
                    if ('---' in infoLine) :
                        rank=0
                        numBars+=1
                    else :
                        rank+=1

                    if rank < 6 and rank > 0:
                        varName=filter(lambda x : x != "",infoLine.split(' '))[1]
                        if varName not in varCounts[curPt]: varCounts[curPt][varName]=0
                        varCounts[curPt][varName] += (6 - rank)

                    if numBars == 2 :
                        rank=0
                        foundTable=False


summaryOut=open(options.output+'/SummaryRankings_%s.txt'%(options.outname),'w')

summaryOut.write("Generated with:\n")
summaryOut.write("anaSubs: %s\n"%options.anasubs)
summaryOut.write("trainings: %s\n"%options.trainings)
summaryOut.write("sigs: %s\n"%options.sigs)
summaryOut.write("bkgs: %s\n"%options.bkgs)
summaryOut.write("trees: %s\n\n\n"%options.trees)

summaryOut.write("p_T 1 TeV:\n")
summaryOut.write("Rank   Variable                        Score\n")
summaryOut.write("-"*50+"\n")
outVarArrpt1=[]
for i in range(1,len(varCounts['pt1'])+1) :
    maxRank=-1
    maxVar =''
    for key in varCounts['pt1'] :
        if maxRank > varCounts['pt1'][key] : continue
        if key in outVarArrpt1: continue

        maxRank = varCounts['pt1'][key]
        maxVar=key

    outVarArrpt1 += [maxVar]

for i in range(0,len(outVarArrpt1)) :
    summaryOut.write('%4s  %30s  %i\n'%('%i'%(i+1),
        outVarArrpt1[i],
        varCounts['pt1'][outVarArrpt1[i]]))
summaryOut.write("-"*50)
summaryOut.write("\n"*5)



summaryOut.write("p_T 5 TeV:\n")
summaryOut.write("Rank   Variable                        Score\n")
summaryOut.write("-"*50+"\n")
outVarArrpt5=[]
for i in range(1,len(varCounts['pt5'])+1) :
    maxRank=-1
    maxVar=''
    for key in varCounts['pt5'] :
        if maxRank > varCounts['pt5'][key] : continue
        if key in outVarArrpt5: continue

        maxRank = varCounts['pt5'][key]
        maxVar=key

    outVarArrpt5 += [maxVar]

for i in range(0,len(outVarArrpt5)) :
    summaryOut.write('%4s  %30s  %i\n'%('%i'%(i+1),
        outVarArrpt5[i],
        varCounts['pt5'][outVarArrpt5[i]]))
summaryOut.write("-"*50)
