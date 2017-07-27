#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import pprint

pp = pprint.PrettyPrinter(indent=2)

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--indir',action="store",type="string",dest="indir",default="./analysis/")
parser.add_option('--outdir',action="store",type="string",dest="outdir",default="./")
parser.add_option('--pts',action="store",type="string",dest="pts",default="pt1,pt5")
parser.add_option('--trainings',action="store",type="string",dest="trainings",default="all,massonly,shapesonly")
parser.add_option('--anasubs',action="store",type="string",dest="anasubs",default="r1_h022_e0175_t220,r05_h01_e005_t220,r05_h002_e001_t220")
parser.add_option('--infos',action="store",type="string",dest="infos",default="tracks,tragam,allpar")
parser.add_option('--procs',action="store",type="string",dest="procs",default="W,Z,g,q,t")

(options, args) = parser.parse_args()

anaSubList   = options.anasubs.split(',')
trainingList = options.trainings.split(',')
ptList       = options.pts.split(',')
procList     = options.procs.split(',')
infoList     = options.infos.split(',')

# prepare arrays
varCounts={}
for anaSub in anaSubList :
    varCounts[anaSub] = {}
    for training in trainingList :
        varCounts[anaSub][training] = {}
        for pt in ptList :
            varCounts[anaSub][training][pt] = {}
            for sig in procList :
                varCounts[anaSub][training][pt][sig] = {}
                for bkg in procList :
                    if sig == bkg : continue
                    varCounts[anaSub][training][pt][sig][bkg] = {}
                    for info in infoList :
                        varCounts[anaSub][training][pt][sig][bkg][info] = {}

if __name__ == '__main__':

    # loop over configurations
    for thisPt,thisAnasub,thisTraining,thisInfo in [(a,b,c,d)
            for a in ptList
            for b in anaSubList
            for c in trainingList
            for d in infoList] :

        # get filenames
        fileRegexp="/trainings_%s_%s/out_*%s*%s*.stdout"%(thisAnasub,thisTraining,thisPt,thisInfo)
        onlyfiles = glob.glob(options.indir+fileRegexp)

        print fileRegexp

        for fname in onlyfiles :
            print '\t - ',fname

            thisSig = (fname.split('out_')[1])[0]
            thisBkg = ((fname.split('out_')[1]).split('_')[1])[0]

            with open(fname,'r') as fIn :
                counter=0
                for line in fIn :
                    if counter > 1 :
                        counter -= 1;
                        continue
                    elif counter == 1 :
                        try :
                            sepVal = ((line.split('|')[1])[1:6])
                            if "0.000" in sepVal :
                                counter == 0
                                continue
                            varCounts[thisAnasub][thisTraining][thisPt][thisSig][thisBkg][thisInfo] = sepVal
                            varCounts[thisAnasub][thisTraining][thisPt][thisBkg][thisSig][thisInfo] = sepVal
                        except IndexError :
                            print "-"*80
                            print "Error with file",fname
                            print "\t\t Problem line: "
                            print "\t\t\t ",line
                            print "-"*80

                    elif "Sepa-" in line :
                        counter = 3
                        continue

# print our data
pp.pprint(varCounts)
print("\n\n")
print("-"*50)
print("Making output file")
print("-"*50)

# open output file
summaryOut=open(options.outdir+'/BDTSeps.tex','w')

# write in tex tables
for thisAnasub in varCounts :
    for thisTraining in varCounts[thisAnasub] :
        summaryOut.write("\n\n\\vspace{2cm}\n")
        summaryOut.write(thisAnasub.replace('_','\_')+', '+thisTraining+'\n')
        summaryOut.write("\\begin{tabular}{|"+('l'*(1 + 1 + len(infoList)))+"|}\n")
        summaryOut.write(" & & " + " & ".join(infoList) + " \\\\\\hline\\hline \n")

        for thisPt in ptList :
            for iSig in range(len(procList)) :
                for iBkg in range(len(procList)) :
                    if iBkg <= iSig : continue

                    thisSig = procList[iSig]
                    thisBkg = procList[iBkg]

                    curLine = ' %s & %s vs. %s '%(thisPt,thisSig,thisBkg)
                    for thisInfo in infoList :
                        curLine = curLine + ' & %s'%varCounts[thisAnasub][thisTraining][thisPt][thisSig][thisBkg][thisInfo]

                    curLine = curLine + ' \\\\\\hline \n'
                    summaryOut.write(curLine)


        summaryOut.write("\\end{tabular}\n")

summaryOut.close()
