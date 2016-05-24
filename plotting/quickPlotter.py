#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import string
import ROOT

import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.12);
ROOT.gStyle.SetPadRightMargin(0.25);
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
def main():
	
	fname1 = "../processing/testdat/test-qq-1k.root"
	fname2 = "../processing/testdat/test-gg-1k.root"

	f1 = ROOT.TFile(fname1);
	f2 = ROOT.TFile(fname2);

	t1_tracks = f1.Get("t_tracks");
	t1_tragam = f1.Get("t_tragam");
	t1_allpar = f1.Get("t_allpar");
	t2_tracks = f2.Get("t_tracks");
	t2_tragam = f2.Get("t_tragam");
	t2_allpar = f2.Get("t_allpar");

	h1_tracks = [];
	h1_tragam = [];
	h1_allpar = [];
	h2_tracks = [];
	h2_tragam = [];
	h2_allpar = [];
	declareHistograms(h1_tracks,"h1_tracks");
	declareHistograms(h1_tragam,"h1_tragam");
	declareHistograms(h1_allpar,"h1_allpar");
	declareHistograms(h2_tracks,"h2_tracks");
	declareHistograms(h2_tragam,"h2_tragam");
	declareHistograms(h2_allpar,"h2_allpar");

	for i in range(t1_tracks.GetEntries()):
		t1_tracks.GetEntry(i);
		t1_tragam.GetEntry(i);
		t1_allpar.GetEntry(i);
		fillHistograms( h1_tracks, t1_tracks );
		fillHistograms( h1_tragam, t1_tragam );
		fillHistograms( h1_allpar, t1_allpar );

	for i in range(t2_tracks.GetEntries()):
		t2_tracks.GetEntry(i);
		t2_tragam.GetEntry(i);
		t2_allpar.GetEntry(i);
		fillHistograms( h2_tracks, t2_tracks );
		fillHistograms( h2_tragam, t2_tragam );
		fillHistograms( h2_allpar, t2_allpar );

	############################################################

	plotsnames = ['jetpt','jetefrac','jeteta',
				  'c1beta0','c1beta1','c1beta2',
				  'multiplicity','mass','mSDbeta0','zlogz'];
	leg1= ['tracks','tracks and photons','all particles'];
	nplots = 10;
	for i in range(nplots):
		makeCanvas( [h1_tracks[i],h1_tragam[i],h1_allpar[i]],
					leg1,
					plotsnames[i]+"-qq-only");
	leg2 = ['q jet','g jet']
	for i in range(nplots):
		makeCanvas( [h1_tracks[i],h2_tracks[i]],
					leg2,
					plotsnames[i]+"-qqvgg-tracks");
	for i in range(nplots):
		makeCanvas( [h1_tragam[i],h2_tragam[i]],
					leg2,
					plotsnames[i]+"-qqvgg-tragam");
	for i in range(nplots):
		makeCanvas( [h1_allpar[i],h2_allpar[i]],
					leg2,
					plotsnames[i]+"-qqvgg-allpar");		

def fillHistograms(hs,t1):
	# print t1_tracks.j_pt[0]
	hs[0].Fill( t1.j_pt[0] );
	hs[1].Fill( t1.j_ptfrac[0] );
	hs[2].Fill( t1.j_eta[0] );
	hs[3].Fill( t1.j_c1_b0[0] );
	hs[4].Fill( t1.j_c1_b1[0] );
	hs[5].Fill( t1.j_c1_b2[0] );
	hs[6].Fill( t1.j_multiplicity[0] );
	hs[7].Fill( t1.j_mass[0] );
	hs[8].Fill( t1.j_mass_mmdt[0] );
	hs[9].Fill( t1.j_zlogz[0] );
	

def declareHistograms(hs,tag):
	hs.append( ROOT.TH1F( "j_pt"+tag, "; jet pT (GeV);",20,0,1500) );
	hs.append( ROOT.TH1F( "j_efrac"+tag, "; jet energy fraction;",20,0,2) );
	hs.append( ROOT.TH1F( "j_eta"+tag, "; eta;",20,-3,3) );
	
	hs.append( ROOT.TH1F( "j_c1_b0"+tag, "; C_{1}^{#beta=0};",20,0,0.5) );
	hs.append( ROOT.TH1F( "j_c1_b1"+tag, "; C_{1}^{#beta=1};",20,0,0.5) );
	hs.append( ROOT.TH1F( "j_c1_b2"+tag, "; C_{1}^{#beta=2};",20,0,0.5) );
	hs.append( ROOT.TH1F( "j_multiplicity"+tag, "; multiplicity;",20,0,200) );

	hs.append( ROOT.TH1F( "j_mass"+tag, "; mass (GeV);",20,0,300) );
	hs.append( ROOT.TH1F( "j_mmdt"+tag, "; m_{SD}^{#beta=0} (GeV);",20,0,300) );
	hs.append( ROOT.TH1F( "j_zlogz"+tag, "; #Sigma z logz;",20,-10,2) );

def makeCanvas(hs,legs,name):
	print "name = ", name;
	colors = [1,2,4,6,7];
	maxval = -999;
	for h in hs: 
		h.Scale(1./h.Integral());
		h.SetLineWidth(2);
		if h.GetMaximum() > maxval: maxval =  h.GetMaximum()
	leg = ROOT.TLegend(0.2,0.7,0.5,0.9)
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.035);
	i = 0;
	for h in hs:
		leg.AddEntry(h,legs[i],"l")
		i+=1;

	c = ROOT.TCanvas("c","c",1000,800);
	hs[0].SetMaximum(maxval*1.5);
	hs[0].Draw("hist");
	i = 0;
	for h in hs: 
		h.SetLineColor(colors[i])
		h.Draw("histsames");
		i+=1;
	leg.Draw();
	c.SaveAs("plots/"+name+".pdf");
	c.SaveAs("plots/"+name+".png");

########################################################################################################################
if __name__ == '__main__':
	main();


