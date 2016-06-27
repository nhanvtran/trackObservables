#! /usr/bin/env python

import os
import glob
import math
from array import array
import sys
import time

import ROOT

def makeCanvas(hists, names, canname, odir, normalize=False,setLogy=False):
	
	#bin = options.ptbin;
	#directory = "figs_bin"+str(bin);
	colors = [2,4,1,6,7];
	
	leg = ROOT.TLegend(0.4,0.7,0.6,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	for i in range(len(names)):
		hists[i].SetLineColor(colors[i]);
		hists[i].SetLineWidth(2);
		if normalize and hists[i].Integral() > 0: hists[i].Scale(1./hists[i].Integral());
		leg.AddEntry(hists[i], names[i], "l");

	max = -999.;
	for hist in hists:
		if max < hist.GetMaximum(): max = hist.GetMaximum();

	can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
	hists[0].SetMaximum( 1.2*max );
	hists[0].SetMinimum( 1e-5 );    
	hists[0].Draw();    
	for i in range(1,len(hists)):
		hists[i].Draw("sames");
	leg.Draw();
	if setLogy: ROOT.gPad.SetLogy();
	can.SaveAs(odir+"/"+canname+".eps");
	can.SaveAs(odir+"/"+canname+".png");
	can.SaveAs(odir+"/"+canname+".pdf");


def makeROCFromHisto(hists,LtoR=True):

	hsig = hists[0];
	hbkg = hists[1];

	nbins = hsig.GetNbinsX();
	binsize = hsig.GetBinWidth(1);
	lowedge = hsig.GetBinLowEdge(1);

	#print "lowedge: ",lowedge

	hsigIntegral = hsig.Integral();
	hbkgIntegral = hbkg.Integral();

	xval = array('d', [])
	yval = array('d', [])
	ctr = 0;
	effBkgPrev = -9999;
	for i in range(1,nbins+1):

			effBkg = 0;
			effSig = 0;

			if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
			else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;

			if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
			else: effSig = hsig.Integral( 1, i )/hsigIntegral;

			#if not effBkg == 0.: print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig;

			xval.append( effSig );
			yval.append( effBkg );

			#effBkgPrev = effBkg;
			ctr = ctr + 1;

	#print nbins, "and ", ctr
	tg = ROOT.TGraph( nbins, xval, yval );
	tg.SetName( "tg"+hsig.GetName() );
	return tg;