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

parser.add_option('--inputs',action="store",type="string",dest="inputs",default="./tmp/")
parser.add_option('--output',action="store",type="string",dest="output",default="./")

(options, args) = parser.parse_args()

varCounts={}
varCounts['pt1']={}
varCounts['pt5']={}

if __name__ == '__main__':

    onlyfiles = [f for f in os.listdir(options.inputs)]

    for fname in onlyfiles :
        foundTable = False;
        numBars = 0;
        if 'pt5' in fname : continue
        if not 'stdout' in fname : continue
        with open('/'.join([options.inputs,fname]),'r') as fin :
            empty=False
            for i, line in enumerate(fin):
                if i == 0 and line == '' : empty=True
                else : break

            if empty : continue

            print fname;
            fout = open(options.output+'/'+'_'.join(fname.replace('out_','').replace('pt1_','').split('_')[0:-1])+'.txt','w')
            rank=0
            for line in fin:
                infoLine=''.join((line+":").split(':')[1:-1])
                if 'Separation' in infoLine :
                    fout.write('pt 1 TeV\n')
                    foundTable=True
                if foundTable :
                    if ('---' in infoLine) :
                        rank=0
                        numBars+=1
                    else :
                        rank+=1
                    if 'Rank' in infoLine : rank = 0

                    if rank < 6 and rank > 0:
                        varName=filter(lambda x : x != "",infoLine.split(' '))[1]
                        if varName not in varCounts['pt1']: varCounts['pt1'][varName]=0
                        varCounts['pt1'][varName] += (6 - rank)

                    fout.write(infoLine)

                    if numBars == 2 :
                        foundTable=False
                        fout.write('\n')

            newname=""
            for tfname in onlyfiles :
                if '_'.join(fname.replace('pt1','pt5').split('_')[0:-1]) in '_'.join(tfname.split('_')[0:-1]) :
                    if 'stdout' in tfname :
                        newname=tfname

            print " -  newname %s "%newname

            foundTable = False;
            numBars = 0;
            with open('/'.join([options.inputs,newname]),'r') as fin5 :
                rank=0
                for line in fin5:
                    infoLine=''.join((line+":").split(':')[1:-1])
                    if 'Separation' in infoLine :
                        fout.write('\n\npt 5 TeV\n')
                        foundTable=True
                    if foundTable :
                        if ('---' in infoLine) :
                            rank=0
                            numBars+=1
                        else :
                            rank+=1
                        if 'Rank' in infoLine : rank = 0

                        if rank < 6 and rank > 0:
                            varName=filter(lambda x : x != "",infoLine.split(' '))[1]
                            if varName not in varCounts['pt5'] : varCounts['pt5'][varName]=0
                            varCounts['pt5'][varName] += (6 - rank)

                        fout.write(infoLine)

                        if numBars == 2 :
                            foundTable=False
                            fout.write('\n')


summaryOut=open(options.output+'/SummaryRankings.txt','w')

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
