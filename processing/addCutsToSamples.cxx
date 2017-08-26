#include <TFile.h>
#include <TObject.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <vector>

void addCutsToSamples(TString fileName, Double_t cutMin, Double_t cutMax) {
    
        // Get allpar information
        TFile* fIn = new TFile(fileName, "UPDATE");
        TTree* t_allpar = ((TTree*) fIn->Get("t_allpar"));
        TTree* t_tragam = ((TTree*) fIn->Get("t_tragam"));
        TTree* t_tracks = ((TTree*) fIn->Get("t_tracks"));

        std::vector<float> * aj_pt=0;
        t_allpar->SetBranchAddress("j_pt", &aj_pt);
        
        Int_t pCutBool = 0;
        Int_t pCutBool_tragam = 0;
        Int_t pCutBool_tracks = 0;
        
        TBranch* passCuts = t_allpar->Branch("j_passCut", &pCutBool, "j_passCut/I");
        TBranch* passCuts_tragam = t_tragam->Branch("j_passCut", &pCutBool_tragam, "j_passCut/I");
        TBranch* passCuts_tracks = t_tracks->Branch("j_passCut", &pCutBool_tracks, "j_passCut/I");
        std::cout<<"Produced new branches."<<std::endl;
        
        for(Int_t i=0; i<t_allpar->GetEntries(); i++) {
            t_allpar->GetEvent(i);
        
            pCutBool = (cutMin <= (*aj_pt)[0]) and ((*aj_pt)[0] <= cutMax);
            std::cout<<"Processed event "<<i<<" with cut pass "<<pCutBool<<std::endl;
        
            pCutBool_tragam=pCutBool;
            pCutBool_tracks=pCutBool;
        
            passCuts->Fill();
        
            t_tragam->GetEvent(i);
            passCuts_tragam->Fill();
        
            t_tracks->GetEvent(i);
            passCuts_tracks->Fill();
        }
        
        t_allpar->Write("t_allpar",TObject::kOverwrite);
        t_tragam->Write("t_tragam",TObject::kOverwrite);
        t_tracks->Write("t_tracks",TObject::kOverwrite);
        
        fIn->Close();
}
