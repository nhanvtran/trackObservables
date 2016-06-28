#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>

TLegend *leg  = new TLegend(0.15625,0.621654,0.4765625,0.803839,NULL,"brNDC");
TCanvas *cROC = new TCanvas("cROC","cROC",700,700);
bool isFirst = true;

// recursion base case: draw and save
void plotROC() {
   leg->Draw();
   cROC->Print("roc_.png");
   cROC->SetLogy();
   cROC->Print("roc_log.png");
   isFirst = true;
}

template<typename... Args>
void plotROC(TString input, TString label, Args... moreLabels) { 

    //TFile * _file0 = new TFile("TMVA_sL_optimized.rooVt");
	//TFile * _file1 = new TFile ("TMVA_QCD_BBvsGSP_fat.root");
    if(isFirst) {
	    // base style
        gROOT->SetStyle("Plain");
	    gStyle->SetPadGridX(0);
	    gStyle->SetPadGridY(0);
	    gStyle->SetOptStat(0);

	    cROC->SetTickx(1);
	    cROC->SetTicky(1);

        // legend
	    leg->SetBorderSize(0);
	    leg->SetTextSize(0.035);
	    leg->SetTextFont(42);
	    leg->SetLineColor(1);
	    leg->SetLineStyle(1);
	    leg->SetLineWidth(1);
	    leg->SetFillColor(0);
	    leg->SetFillStyle(0);

        // plot first ROC
	    TFile * file = new TFile(input);
	    TH1D* MVA_BDTG_effBvsS = (TH1D *) file->Get("Method_BDT/BDTG/MVA_BDTG_effBvsS");
	    MVA_BDTG_effBvsS->SetTitle("");
	    MVA_BDTG_effBvsS->SetLineColor(kBlue+1);
        MVA_BDTG_effBvsS->SetLineWidth(3);
	    
	    MVA_BDTG_effBvsS->GetXaxis()->SetTitle("Signal efficiency");
        MVA_BDTG_effBvsS->GetYaxis()->SetTitle("Background efficiency");
	    MVA_BDTG_effBvsS->GetXaxis()->SetTitleSize(0.04);
	    MVA_BDTG_effBvsS->GetYaxis()->SetTitleSize(0.04);
	    MVA_BDTG_effBvsS->GetXaxis()->SetTitleOffset(1.05);
	    MVA_BDTG_effBvsS->GetYaxis()->SetTitleOffset(1.2);
	    MVA_BDTG_effBvsS->GetXaxis()->SetLabelSize(0.03);
	    MVA_BDTG_effBvsS->GetYaxis()->SetLabelSize(0.03);
        
        MVA_BDTG_effBvsS->Draw();
	    
        leg->AddEntry(MVA_BDTG_effBvsS, label,"l");
    } else { 
		TFile * file = new TFile(input);
	    TH1D* MVA_BDTG_effBvsS = (TH1D *) file->Get("Method_BDT/BDTG/MVA_BDTG_effBvsS");
		MVA_BDTG_effBvsS->SetLineColor(kGreen+2);
		MVA_BDTG_effBvsS->SetLineWidth(3);
		MVA_BDTG_effBvsS->Draw("same");
		leg->AddEntry(MVA_BDTG_effBvsS, label, "l");
	}

    // continue to call while there are extra args
    isFirst = false;
    plotROC(moreLabels...);

}
