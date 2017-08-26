#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from TMVAhelper import *
from utilities import *

# ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
# ROOT.setTDRStyle();
# ROOT.gStyle.SetPadTopMargin(0.06);
# ROOT.gStyle.SetPadLeftMargin(0.16);
# ROOT.gStyle.SetPadRightMargin(0.10);
# ROOT.gStyle.SetPalette(1);
# ROOT.gStyle.SetPaintTextFormat("1.1f");


class observableContainer:

	# -------------------------------------
	def __init__(self,signame,bkgname,variables,cuts,label,treeName,weightloc='weights',MVAMethod="BDTG"):

		self._f_sig = ROOT.TFile(signame);
		self._f_bkg = ROOT.TFile(bkgname);
		self._t_sig = self._f_sig.Get(treeName);
		self._t_bkg = self._f_bkg.Get(treeName);

		self._discVariables = variables;#,"MHT","nJets_30"];

		#tmva cut
		self._cutstring = "(";
		cutctr = 0;
		for cut in cuts:
			self._cutstring += "(" + cut[0] + " >= " + str(cut[1]) + ") && (" + cut[0] + " <= " + str(cut[2]) + ")"
			if cutctr < len(cuts) - 1: self._cutstring += "&&";
			cutctr+=1;
		self._cutstring += ")";
		print "cutstring = ", self._cutstring;

		# self._spectatorVariables = ["lheWeight"];
		self._spectatorVariables = [];
		self._trees = [self._t_sig,self._t_bkg];
		self._bdt = TMVAhelper("MVA_"+label,self._discVariables,self._trees,self._spectatorVariables,weightloc);
		print ";;weightloc = ",weightloc
	# -------------------------------------
	def doTraining(self):
		self._bdt.train(self._cutstring);

	def readMVA(self,method="BDTG"):
		# read in weight files (one can train and loop at different times)
		self._bdt.read(method);

	def evaluateMVA(self,val,method="BDTG"):
		return self._bdt.evaluate(val,method);



