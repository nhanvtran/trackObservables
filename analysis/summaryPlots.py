#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from observableContainer import *
from TMVAhelper import *
from utilities import *

ROOT.gROOT.ProcessLine(".L tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadTopMargin(0.06);
ROOT.gStyle.SetPadLeftMargin(0.06);
ROOT.gStyle.SetPadRightMargin(0.06);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

(options, args) = parser.parse_args()

############################################################

# observableTraining takes in 2 root files, list of observables, spectator observables ... launches a CONDOR job
# TMVAhelper.py is used by observableTraining
# analysis.py defines the list of trainings, the signal and background process

########################################################################################################################
########################################################################################################################
def getFilesRecursively(dir,searchstring,additionalstring = None):
	
	cfiles = [];
	for root, dirs, files in os.walk(dir):
		for file in files:	
			if searchstring in file:
				if additionalstring == None or additionalstring in file:
					cfiles.append(os.path.join(root, file))
	return cfiles;

def fileVariables(fname):
    variables = [];
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
            variables.append( l.strip().split()[0] );
    #return i + 1	
    return variables;

def findEff(name,fn,effval):
	index = -1;
	if effval == 10: index = 1;
	if effval == 25: index = 2;

	thevalue = 0;
	curf = open(fn,'r');
	for line in curf:
		curline = line.strip().split();
		if curline[0] == name: thevalue = curline[index];
	return thevalue;


def makeSummary( mass, background, effval):

	plotDir   = './training/plots';
	#plotDir   = '../plots';
	searchtag = "%s__%s" % (str(mass),background);
	filesToRead = getFilesRecursively(plotDir,"RocSummary",searchtag);
	print filesToRead, len(filesToRead); 

	signals = [];
	masses = [str(mass)];
	for i in range(1,5): 
		tag1,tag2 = "G","G";
		for ai in range(1,i+1): tag1+="j";
		for ai in range(1,i+1): tag2+="j";
		for m in masses: 
			signals.append( tag1+"N1_"+tag2+"N1__"+m+"_" );
			if i < 4: signals.append( tag1+"N1_"+tag2+"jN1__"+m+"_" );

	
	variables = fileVariables(filesToRead[0]);
	nvariables = len(variables);
	print signals, nvariables
	histograms = [];
	for i in range(nvariables):
		histograms.append(ROOT.TH1F("h_"+variables[i],";n partons;1/#epsilon_{bkg}",len(signals),2,len(signals)+2));

	theMax = -99;
	for j in range(len(signals)): 
		for fn in filesToRead:
			if signals[j] in fn: 
				print signals[j], "------"
				for i in range(nvariables):
					curvalue = float(findEff(variables[i],fn,effval));
					print variables[i],curvalue,j
					histograms[i].SetBinContent(j+1,curvalue);
					if curvalue > theMax: theMax = curvalue;

	cSum = ROOT.TCanvas("cSum","cSum",1000,800);
	
	leg = ROOT.TLegend(0.2,0.75,0.8,0.9);
	leg.SetNColumns(2);
	leg.SetFillColor(0);
	leg.SetBorderSize(0);

	colors = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,12,12]
	style =  [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2, 1, 2]
	histograms[0].SetMaximum(theMax*1.25);
	histograms[0].SetMinimum(0.);
	for i in range(len(histograms)):
		histograms[i].SetLineColor(colors[i]);
		histograms[i].SetLineWidth(4);
		histograms[i].SetLineStyle(style[i]);

		if i == 0: histograms[i].Draw("c");
		else: histograms[i].Draw("csames");
		leg.AddEntry(histograms[i],variables[i],"l");
	leg.Draw();

	banner = ROOT.TLatex(0.18,0.96,("m_{Gluino} = "+str(mass)+", background = "+background+", #epsilon_{sig} = "+str(effval)+"%"));
	banner.SetNDC()
	banner.SetTextSize(0.04)
	banner.Draw();
	tag = "%s_%s_%s" % (str(mass),background,str(effval));
	cSum.SaveAs("summaryplots/summary_"+tag+".png");
	cSum.SaveAs("summaryplots/summary_"+tag+".pdf");
	histograms[0].SetMaximum(theMax*10.);
	histograms[0].SetMinimum(2.);
	ROOT.gPad.SetLogy();
	cSum.SaveAs("summaryplots/summary_"+tag+"log.png");
	cSum.SaveAs("summaryplots/summary_"+tag+"log.pdf");

#------------------------------------------------#------------------------------------------------
#------------------------------------------------#------------------------------------------------
#------------------------------------------------#------------------------------------------------
#------------------------------------------------#------------------------------------------------

if __name__ == '__main__':

	print "Make summary plots..."

	masses = [500,1000,1500];
	# bkgs   = ['ttbar','QCD','Wjets','znunu'];
	bkgs   = ['ttbar','Wjets','znunu'];

	for mass in masses:
		for bkg in bkgs:
			print "Make mass = ",str(mass),", bkg = ", bkg;
			
			makeSummary( mass, bkg, 10 );
			makeSummary( mass, bkg, 25 );

