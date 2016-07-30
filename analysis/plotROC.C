#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>

// drawing utils
TLegend *leg  = new TLegend(0.15625,0.621654,0.4765625,0.803839,NULL,"brNDC");
TCanvas *cROC = new TCanvas("cROC","cROC",700,700);
bool isFirst = true;

// style utils
std::vector<int> colors = { kBlue+1, kGreen+2, kRed+2, kOrange, kMagenta, kYellow, kViolet, kCyan };
int iColor = 0;
int sColor = 0;
int iCMax = colors.size();
bool dotted = false;

// recursion base case: draw and save
void plotROC() {
    leg->Draw();
    cROC->Print("roc_.png");
    cROC->Print("roc_.pdf");
    cROC->SetLogy();
    cROC->Print("roc_log.png");
    cROC->Print("roc_log.pdf");
    
    leg->Clear();
    cROC->SetLogy(0);
    cROC->Clear();

    iColor = 0;
    sColor = 0;
    dotted = false;

    isFirst = true;
}

void plotROC(TString filename) {
    leg->Draw();
    cROC->Print("roc_"+filename+".png");
    cROC->Print("roc_"+filename+".pdf");
    cROC->SetLogy();
    cROC->Print("roc_log_"+filename+".png");
    cROC->Print("roc_log_"+filename+".pdf");
    
    leg->Clear();
    cROC->SetLogy(0);
    cROC->Clear();

    iColor = 0;
    sColor = 0;
    dotted = false;

    isFirst = true;
}

template<typename... Args>
void plotROC(TString input, TString label, Args... moreLabels) { 

    dotted = ( sColor % 2 != 0 );

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
        if(!((bool) file->GetListOfKeys()->At(0))) goto justPlotNext;
        TH1F* MVA_BDTG_effBvsS = (TH1F *) file->Get("Method_BDT/BDTG/MVA_BDTG_effBvsS");
        MVA_BDTG_effBvsS->SetTitle("");
        MVA_BDTG_effBvsS->SetLineColor(colors.at(0));
        MVA_BDTG_effBvsS->SetLineStyle(kSolid); 
        MVA_BDTG_effBvsS->SetLineWidth(3);

        MVA_BDTG_effBvsS->GetXaxis()->SetTitle("Signal efficiency");
        MVA_BDTG_effBvsS->GetYaxis()->SetTitle("Background efficiency");
        MVA_BDTG_effBvsS->GetXaxis()->SetTitleSize(0.04);
        MVA_BDTG_effBvsS->GetYaxis()->SetTitleSize(0.04);
        MVA_BDTG_effBvsS->GetXaxis()->SetTitleOffset(1.05);
        MVA_BDTG_effBvsS->GetYaxis()->SetTitleOffset(1.2);
        MVA_BDTG_effBvsS->GetXaxis()->SetLabelSize(0.03);
        MVA_BDTG_effBvsS->GetYaxis()->SetLabelSize(0.03);
        MVA_BDTG_effBvsS->SetMaximum(1);

        MVA_BDTG_effBvsS->Draw();

        leg->AddEntry(MVA_BDTG_effBvsS, label,"l");
        //file->Close();
    } else { 
        TFile * file = new TFile(input);
        if(!((bool) file->GetListOfKeys()->At(0))) goto justPlotNext;
        TH1F* MVA_BDTG_effBvsS = (TH1F *) file->Get("Method_BDT/BDTG/MVA_BDTG_effBvsS");
        MVA_BDTG_effBvsS->SetLineColor(colors.at(iColor)+0*sColor);
        if(dotted) { 
            MVA_BDTG_effBvsS->SetLineStyle(kDotted); 
            MVA_BDTG_effBvsS->SetMarkerStyle(kDot); 
            MVA_BDTG_effBvsS->SetMarkerSize(4); 
        } else { 
            MVA_BDTG_effBvsS->SetLineStyle(kSolid); 
        }
        MVA_BDTG_effBvsS->SetLineWidth(3);
        MVA_BDTG_effBvsS->SetMaximum(1);
        MVA_BDTG_effBvsS->Draw("same");
        leg->AddEntry(MVA_BDTG_effBvsS, label, (dotted ? "lp" : "l"));
        //file->Close();
    }

    iColor++;

    if(iColor == iCMax) {
        iColor=0;
        sColor += 1;
    } 

    // continue to call while there are extra args
    isFirst = false;

    justPlotNext:
        plotROC(moreLabels...);

}
