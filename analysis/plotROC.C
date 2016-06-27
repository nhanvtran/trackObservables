#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>


void plotROC(TString input1, TString label1, TString input2 , TString label2,  TString input3,  TString label3, TString input4, TString label4,  TString input5,  TString label5 , TString input6,  TString label6 ){


	gROOT->SetStyle("Plain");
	gStyle->SetPadGridX(0);
	gStyle->SetPadGridY(0);
	gStyle->SetOptStat(0);
	TCanvas *cROC = new TCanvas("cROC","cROC", 700, 700);

	cROC->SetTickx(1);
	cROC->SetTicky(1);
	TLegend *leg = new TLegend(0.15625,0.621654,0.4765625,0.803839,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);




	//TFile * _file0 = new TFile("TMVA_sL_optimized.root");
	//TFile * _file1 = new TFile ("TMVA_QCD_BBvsGSP_fat.root");
	TFile * file1 = new TFile(input1);
	file1->cd();
	TDirectoryFile * d = (TDirectoryFile *)file1->Get("Method_BDT/BDTG");
	d->cd();
	//cROC->cd();
	MVA_BDTG_effBvsS->SetTitle("");
	MVA_BDTG_effBvsS->SetLineColor(kBlue+1);
        MVA_BDTG_effBvsS->SetLineWidth(3);
	MVA_BDTG_effBvsS->GetYaxis()->SetTitle("Background efficiency");
	MVA_BDTG_effBvsS->GetYaxis()->SetTitleSize(0.04);
	MVA_BDTG_effBvsS->GetXaxis()->SetTitle("Signal efficiency");
	MVA_BDTG_effBvsS->GetXaxis()->SetTitleSize(0.04);
	MVA_BDTG_effBvsS->GetXaxis()->SetTitleOffset(1.05);
	MVA_BDTG_effBvsS->GetYaxis()->SetTitleOffset(1.2);
	MVA_BDTG_effBvsS->GetXaxis()->SetLabelSize(0.03);
	MVA_BDTG_effBvsS->GetYaxis()->SetLabelSize(0.03);
        MVA_BDTG_effBvsS->Draw();
	leg->AddEntry(MVA_BDTG_effBvsS, label1,"l");
	if(input2!=""){
		TFile * file2 = new TFile(input2);
		file2->cd();
		TDirectoryFile * d2 = file2->Get("Method_BDT/BDTG");
		d2->cd();
		MVA_BDTG_effBvsS->SetLineColor(kGreen+2);
		MVA_BDTG_effBvsS->SetLineWidth(3);
		MVA_BDTG_effBvsS->Draw("same");
		leg->AddEntry(MVA_BDTG_effBvsS, label2, "l");

	}
	if(input3!=""){
		TFile * file3 = new TFile(input3);
		file3->cd();
		TDirectoryFile * d3 = file3->Get("Method_BDT/BDTG");
		d3->cd();
		MVA_BDTG_effBvsS->SetLineColor(kRed+2);
		MVA_BDTG_effBvsS->SetLineWidth(3);
		MVA_BDTG_effBvsS->Draw("same");
		leg->AddEntry(MVA_BDTG_effBvsS, label3,"l");
	
	}
	if(input4!=""){
		TFile * file4 = new TFile(input4);
		file4->cd();
		TDirectoryFile * d4 = file4->Get("Method_BDT/BDTG");
		d4->cd();
		MVA_BDTG_effBvsS->SetLineColor(kBlue+1);
		MVA_BDTG_effBvsS->SetLineWidth(3);
		MVA_BDTG_effBvsS->SetLineStyle(2);
		MVA_BDTG_effBvsS->Draw("same");
		leg->AddEntry(MVA_BDTG_effBvsS, label4);
	
	}
	if(input5!=""){
		TFile * file5 = new TFile(input5);
		file5->cd();
		TDirectoryFile * d5 = file5->Get("Method_BDT/BDTG");
		d5->cd();
		 MVA_BDTG_effBvsS->SetLineColor(kGreen+2);
                MVA_BDTG_effBvsS->SetLineWidth(3);
                MVA_BDTG_effBvsS->SetLineStyle(2);
		MVA_BDTG_effBvsS->Draw("same");
		leg->AddEntry(MVA_BDTG_effBvsS, label5);

	}
	if(input6!=""){
                TFile * file6 = new TFile(input6);
                file6->cd();
                TDirectoryFile * d6 = file6->Get("Method_BDT/BDTG");
                d6->cd();
		MVA_BDTG_effBvsS->SetLineColor(kRed+2);
                MVA_BDTG_effBvsS->SetLineWidth(3);
                MVA_BDTG_effBvsS->SetLineStyle(2);

                MVA_BDTG_effBvsS->Draw("same");
                leg->AddEntry(MVA_BDTG_effBvsS, label6);

        }
	leg->Draw();
	cROC->Print("roc_.png");
	cROC->SetLogy();
	cROC->Print("roc_log.png");
}

