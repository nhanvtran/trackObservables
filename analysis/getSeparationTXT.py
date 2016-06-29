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
            fout = open('_'.join(fname.replace('out_','').replace('pt1_','').split('_')[0:-1])+'.txt','w')
            for line in fin:
                infoLine=''.join((line+":").split(':')[1:-1])
                if 'Separation' in infoLine :
                    fout.write('pt 1 TeV\n')
                    foundTable=True
                if foundTable :
                    if ('---' in infoLine) : numBars+=1

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
                for line in fin5:
                    infoLine=''.join((line+":").split(':')[1:-1])
                    if 'Separation' in infoLine :
                        fout.write('\n\npt 5 TeV\n')
                        foundTable=True
                    if foundTable :
                        if ('---' in infoLine) : numBars+=1

                        fout.write(infoLine)

                        if numBars == 2 :
                            foundTable=False
                            fout.write('\n')


