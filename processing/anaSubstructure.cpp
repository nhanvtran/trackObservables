#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TMath.h"
#include <cmath>
#include "TTree.h"
#include <iostream>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/DistanceMeasure.hh"


#include "LHEF.h"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<float>+;
//#endif

// using namespace std;
// using namespace fastjet;
// using namespace fastjet::contrib;

/*
 
 TO COMPILE:
 
 export ROOTSYS=~/Desktop/root
 export PATH=$ROOTSYS/bin:$PATH
 export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH 
 
 c++ -o anaSubstructure `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaSubstructure.cpp
 
 TO RUN:     
 
 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 
 
 */

//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
int evtCtr;
float RPARAM;

int njets;
std::vector<float> j_ptfrac;
std::vector<float> j_pt;
std::vector<float> j_eta;
std::vector<float> j_mass;
std::vector<float> j_tau1_b1;
std::vector<float> j_tau2_b1;
std::vector<float> j_tau3_b1;
std::vector<float> j_tau1_b2;
std::vector<float> j_tau2_b2;
std::vector<float> j_tau3_b2;
std::vector<float> j_tau21_b2;
std::vector<float> j_tau21_b1;
std::vector<float> j_tau32_b2;
std::vector<float> j_tau32_b1;
std::vector<float> j_zlogz;
std::vector<float> j_c1_b0;
std::vector<float> j_c1_b1;
std::vector<float> j_c1_b2;
std::vector<float> j_c2_b1;
std::vector<float> j_c2_b2;
std::vector<float> j_d2_b1;
std::vector<float> j_d2_b2;
std::vector<float> j_qjetVol;
std::vector<float> j_mass_trim;
std::vector<float> j_mass_mmdt;
std::vector<float> j_mass_prun;
std::vector<float> j_mass_sdb2;
std::vector<float> j_mass_sdm1;
std::vector<float> j_multiplicity;

// 211, -211: PI+, PI-
// 3122, -3122: Lambda 
// 22: photons
// 130: K0L
// 310: K0S
// 321, -321: K+,K-
// 2212: proton
// 2112: neutron
// 3112, 3222: Sigma-, Sigma+
// 3322: Xi0
// 3312: Xi-     

std::map< TString,std::vector<int> > ids {
    { "c" , {211,-211,321,-321,2212,-2212,3112,-3112,3222,-3222,3312,-3312} },
    { "p" , {22,111}                                                        },
    { "n" , {3122,-3122,130,310,2112,-2112,3322,-3322}                      },
    { "l" , {-11,11,-13,13,-15,15}                                          }
};

std::vector<int> c_ids = ids["c"];
std::vector<int> p_ids = ids["p"]; 
std::vector<int> n_ids = ids["n"]; 
std::vector<int> l_ids = ids["l"]; 

int activeAreaRepeats = 1;
double ghostArea = 0.01;
double ghostEtaMax = 7.0;
    
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, TTree* t_tracks, TTree* t_tragam, TTree* t_allpar);
void declareBranches(TTree* t);
void clearVectors();
void PushBackJetInformation(fastjet::PseudoJet jet, int particleContentFlag);

// smearing functions
void smearJetPt(fastjet::PseudoJet jet);
std::vector<fastjet::PseudoJet> discretizeJet(fastjet::PseudoJet jet, bool clusterJet = false);

////////////////////-----------------------------------------------

int main (int argc, char **argv) {
    
    std::string type  = argv[1];   // type "gg" or "qq"
    std::string indir = argv[2];   // where to find input files 
    int min = atoi(argv[3]);      // events to run over
    int max = atoi(argv[4]);      // events to run over
    int tag = atoi(argv[5]);      // index for condorizing 
    
    RPARAM = 0.8;

    char inName[192];
    sprintf( inName, "%s/%s.lhe",indir.c_str(),type.c_str() );
    std::cout << "fname = " << inName << std::endl;
    std::ifstream ifsbkg (inName) ;
    LHEF::Reader reader(ifsbkg) ;

    char outName[192];
    sprintf(outName, "processed-%s-%i.root", type.c_str(), tag);
    TFile *f = TFile::Open(outName,"RECREATE");
    TTree *t_tracks = new TTree("t_tracks","Tree with vectors");
    TTree *t_tragam = new TTree("t_tragam","Tree with vectors");
    TTree *t_allpar = new TTree("t_allpar","Tree with vectors");
    declareBranches(t_tracks);
    declareBranches(t_tragam);
    declareBranches(t_allpar);

    // evtCtr = 0;
    std::vector < fastjet::PseudoJet > particles;

    // loop over events
    while ( reader.readEvent () ) {
        
        ++evtCtr;
        if (evtCtr < min) continue;
        if (evtCtr > max) break;
        
        if (evtCtr % 100 == 0) std::cout << "event " << evtCtr << "\n";
        
        // per event
        particles.clear();

        // std::cout << "reader.hepeup.IDUP.size() = " << reader.hepeup.IDUP.size() << std::endl;
        for (unsigned int i = 0 ; i < reader.hepeup.IDUP.size(); ++i){

            if (reader.hepeup.ISTUP.at(i) == 1){
                float px = reader.hepeup.PUP.at(i).at(0);
                float py = reader.hepeup.PUP.at(i).at(1);
                float pz = reader.hepeup.PUP.at(i).at(2);
                float e  = reader.hepeup.PUP.at(i).at(3);                                    
                fastjet::PseudoJet curpar   = fastjet::PseudoJet( px, py, pz, e );
                int pdgid = reader.hepeup.IDUP.at(i);
                curpar.set_user_index( pdgid );
                particles.push_back( curpar );
            }   

        }

        analyzeEvent( particles, t_tracks, t_tragam, t_allpar );        
    }
    
    std::cout << "finish loop" << std::endl;
    
    f->cd();
    t_tracks->Write();
    t_tragam->Write();
    t_allpar->Write();    
    f->Close();
    
    delete f;
    return 0 ;
}

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, TTree* t_tracks, TTree* t_tragam, TTree* t_allpar){
    
    // recluster on the fly....
    fastjet::JetDefinition   jetDef(fastjet::antikt_algorithm, RPARAM);    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(100.0));

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    njets = out_jets.size();
    // std::cout << "Number of jets in the event = " << njets << std::endl;

    // All Particles Looop
    for (unsigned int i = 0; i < out_jets.size(); i++){
        PushBackJetInformation(out_jets[i],0);
    }
    t_allpar->Fill();
    clearVectors();

    for (unsigned int i = 0; i < out_jets.size(); i++){
        PushBackJetInformation(out_jets[i],1);
    }
    t_tracks->Fill();
    clearVectors();

    // std::cout << "tracks" << std::endl;

    for (unsigned int i = 0; i < out_jets.size(); i++){
        PushBackJetInformation(out_jets[i],2);
    }
    t_tragam->Fill();
    clearVectors();
    
    // std::cout << "tracks+photons" << std::endl;

    thisClustering->delete_self_when_unused();
    
    // std::cout << "delete clustering" << std::endl;

}

// ----------------------------------------------------------------------------------
void declareBranches( TTree* t ){

    t->Branch("njets"            , &njets            );
    t->Branch("j_ptfrac"         , &j_ptfrac         );
    t->Branch("j_pt"             , &j_pt             );
    t->Branch("j_eta"            , &j_eta            );
    t->Branch("j_mass"           , &j_mass           );
    t->Branch("j_tau1_b1"        , &j_tau1_b1        );
    t->Branch("j_tau2_b1"        , &j_tau2_b1        );
    t->Branch("j_tau3_b1"        , &j_tau3_b1        );    
    t->Branch("j_tau1_b2"        , &j_tau1_b2        );
    t->Branch("j_tau2_b2"        , &j_tau2_b2        );
    t->Branch("j_tau3_b2"        , &j_tau3_b2        );    
    t->Branch("j_tau21_b1"       , &j_tau21_b1       );
    t->Branch("j_tau21_b2"       , &j_tau21_b2       );
    t->Branch("j_tau21_b1"       , &j_tau21_b1       );
    t->Branch("j_tau21_b2"       , &j_tau21_b2       );
    t->Branch("j_tau32_b1"       , &j_tau32_b1       );
    t->Branch("j_tau32_b2"       , &j_tau32_b2       );
    t->Branch("j_zlogz"          , &j_zlogz          );
    t->Branch("j_c1_b0"          , &j_c1_b0          );
    t->Branch("j_c1_b1"          , &j_c1_b1          );
    t->Branch("j_c1_b2"          , &j_c1_b2          );
    t->Branch("j_c2_b1"          , &j_c2_b1          );
    t->Branch("j_c2_b2"          , &j_c2_b2          );
    t->Branch("j_d2_b1"          , &j_d2_b1          );
    t->Branch("j_d2_b2"          , &j_d2_b2          );
    t->Branch("j_mass_trim"      , &j_mass_trim      );
    t->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    t->Branch("j_mass_prun"      , &j_mass_prun      );
    t->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    t->Branch("j_mass_sdm1"      , &j_mass_sdm1      );
    t->Branch("j_multiplicity"   , &j_multiplicity   );

}

// ----------------------------------------------------------------------------------
void clearVectors(){
    j_pt.clear();
    j_ptfrac.clear();
    j_eta.clear();
    j_mass.clear();
    j_tau1_b1.clear();
    j_tau2_b1.clear();
    j_tau3_b1.clear();    
    j_tau1_b2.clear();
    j_tau2_b2.clear();
    j_tau3_b2.clear();    
    j_tau21_b1.clear();
    j_tau21_b2.clear();
    j_tau32_b1.clear();
    j_tau32_b2.clear();
    j_zlogz.clear();
    j_c1_b0.clear();
    j_c1_b1.clear();
    j_c1_b2.clear();
    j_c2_b1.clear();
    j_c2_b2.clear();
    j_d2_b1.clear();
    j_d2_b2.clear();
    j_qjetVol.clear();
    j_mass_trim.clear();
    j_mass_mmdt.clear();
    j_mass_prun.clear();
    j_mass_sdb2.clear();
    j_mass_sdm1.clear();
    j_multiplicity.clear();
}

// ----------------------------------------------------------------------------------
void PushBackJetInformation(fastjet::PseudoJet jet, int particleContentFlag){


    // recluster on the fly....
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, RPARAM);    
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );

    // std::cout << "old jet pt = " << jet.pt() << ", and particle content flag = " << particleContentFlag << std::endl;

    fastjet::PseudoJet curjet;  
    std::vector< fastjet::PseudoJet > newparticles;
    bool jet_has_no_particles = false;
    if (particleContentFlag == 0) { 
        curjet = discretizeJet(jet,true)[0];
    }
    else if (particleContentFlag > 0) {

        std::vector< fastjet::PseudoJet > discretePar;
        discretePar = discretizeJet(jet,false);
        
        // make a new set of particles
        for (int j = 0; j < discretePar.size(); j++){
            // std::cout << "pdg ids = " << jet.constituents().at(j).user_index() << std::endl;
            if (std::find(c_ids.begin(), c_ids.end(), discretePar.at(j).user_index()) != c_ids.end()) {
                newparticles.push_back(discretePar.at(j));
            }

            if (std::find(p_ids.begin(), p_ids.end(), discretePar.at(j).user_index()) != p_ids.end() 
                    && particleContentFlag == 2) {
                newparticles.push_back(discretePar.at(j));
            }
        }

        // std::cout << "now cluster these particles " << newparticles.size() << std::endl;
        if (newparticles.size() > 0){
            // cluster them
            fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newparticles, jetDef, fjAreaDefinition);
            std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));        
            
            // std::cout << "now fill " << out_jets.size() <<  std::endl;
            // fill into curjet
            curjet = out_jets[0];
        }
        else {
            curjet = fastjet::PseudoJet(0,0,0,0);
            jet_has_no_particles = true;
        }
    }

    // groomers/taggers
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    
    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0;
    fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdb2(2.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdm1(-1.0, zcut_sd, mu_sd);
    
    // n-subjettiness    
    double beta = 1;      // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
    double R0   = RPARAM; // Characteristic jet radius for normalization              
    double Rcut = RPARAM; // maximum R particles can be from axis to be included in jet   
    // beta = 1                   
    fastjet::contrib::Nsubjettiness nSub1KT_b1(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness nSub2KT_b1(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness nSub3KT_b1(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));

    // beta = 2
    fastjet::contrib::Nsubjettiness nSub1KT_b2(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
    fastjet::contrib::Nsubjettiness nSub2KT_b2(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
    fastjet::contrib::Nsubjettiness nSub3KT_b2(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));

    // ECF
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b1(2,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b2(2,2,fastjet::contrib::EnergyCorrelator::pt_R);

    fastjet::contrib::EnergyCorrelator ECF_E1_b1 (1,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelator ECF_E1_b2 (1,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelator ECF_E2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelator ECF_E2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelator ECF_E3_b1 (3,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelator ECF_E3_b2 (3,2.0,fastjet::contrib::EnergyCorrelator::pt_R);

    j_pt.push_back( curjet.pt() );
    j_eta.push_back( curjet.eta() );
    j_mass.push_back( curjet.m() );                
    j_ptfrac.push_back( curjet.pt()/jet.pt() );

    if (!jet_has_no_particles){

        // N-subjettiness
        j_tau1_b1.push_back(  nSub1KT_b1(curjet) );        
        j_tau2_b1.push_back(  nSub2KT_b1(curjet) );        
        j_tau3_b1.push_back(  nSub3KT_b1(curjet) );        
        j_tau1_b2.push_back(  nSub1KT_b2(curjet) );        
        j_tau2_b2.push_back(  nSub2KT_b2(curjet) );  
        j_tau3_b2.push_back(  nSub3KT_b2(curjet) );  
        j_tau21_b1.push_back( nSub2KT_b1(curjet) / nSub1KT_b1(curjet) );
        j_tau21_b2.push_back( nSub2KT_b2(curjet) / nSub1KT_b2(curjet) );
        j_tau32_b1.push_back( nSub3KT_b1(curjet) / nSub2KT_b1(curjet) );
        j_tau32_b2.push_back( nSub3KT_b2(curjet) / nSub2KT_b2(curjet) );
        
        // Z log Z
        float zlogz = 0.;
        for (int i = 0; i < curjet.constituents().size(); i++){
            float conste = curjet.constituents().at(i).e()/curjet.e();
            zlogz += conste * log(conste);
        }

        // energy correlator     
        j_zlogz.push_back( zlogz );   
        j_c1_b0.push_back( C1_b0(curjet) );
        j_c1_b1.push_back( C1_b1(curjet) );
        j_c1_b2.push_back( C1_b2(curjet) );
        j_c2_b1.push_back( C2_b1(curjet) );
        j_c2_b2.push_back( C2_b2(curjet) );

        double cur_e1_b1 = ECF_E1_b1(curjet);
        double cur_e1_b2 = ECF_E1_b2(curjet);
        double cur_e2_b1 = ECF_E2_b1(curjet);
        double cur_e2_b2 = ECF_E2_b2(curjet);
        double cur_e3_b1 = ECF_E3_b1(curjet);
        double cur_e3_b2 = ECF_E3_b2(curjet);
        j_d2_b1.push_back( cur_e3_b1 * pow(cur_e1_b1,3) / pow(cur_e2_b1,3) );
        j_d2_b2.push_back( cur_e3_b2 * pow(cur_e1_b1,3) / pow(cur_e2_b2,3) );
        
        j_mass_trim.push_back( trimmer1( curjet ).m() );
        j_mass_prun.push_back( pruner1( curjet ).m() );    
        j_mass_mmdt.push_back( soft_drop_mmdt( curjet ).m() );
        j_mass_sdb2.push_back( soft_drop_sdb2( curjet ).m() );
        j_mass_sdm1.push_back( soft_drop_sdm1( curjet ).m() );
        
        j_multiplicity.push_back( (float) curjet.constituents().size() );

    }
    else {

        // N-subjettiness
        j_tau1_b1.push_back( -99 );        
        j_tau2_b1.push_back( -99 );        
        j_tau1_b2.push_back( -99 );        
        j_tau2_b2.push_back( -99 );  
        j_tau21_b1.push_back( -99 );
        j_tau21_b2.push_back( -99 );
        
        // energy correlator     
        j_zlogz.push_back( -99 );
        j_c1_b0.push_back( -99 );
        j_c1_b1.push_back( -99 );
        j_c1_b2.push_back( -99 );
        j_c2_b1.push_back( -99 );
        j_c2_b2.push_back( -99 );
        j_d2_b1.push_back( -99 );
        j_d2_b2.push_back( -99 );
        
        j_mass_trim.push_back( -99 );
        j_mass_prun.push_back( -99 );
        j_mass_mmdt.push_back( -99 );
        j_mass_sdb2.push_back( -99 );
        j_mass_sdm1.push_back( -99 );
        
        j_multiplicity.push_back( -99 );

    }
}

// ----------------------------------------------------------------------------------
void smearJetPt(fastjet::PseudoJet jet) {
    static TF1 *smearDist = new TF1("smearGaussian","gaus",0,10000);

    Double_t energyResolution=0;
    for(auto& it: ids) { 
        std::vector<int> tidSet = it.second;

        if (std::find(tidSet.begin(), tidSet.end(), jet.user_index()) != tidSet.end()) {
            if(it.first=="c")                                     // -------------- charged hadron resolution
                // NOTE: minimize neutral hadron vs. track resolution
                energyResolution=sqrt(TMath::Min(pow(0.00001*jet.pt()  ,2)+pow(0.005,2),
                                                 pow(0.38/sqrt(jet.E()),2)+pow(0.01 ,2)));
            else if (it.first=="p")                               // ---------------------- photon resolution
                energyResolution=sqrt(pow(0.10/sqrt(jet.E()),2)+pow(0.0075,2));
            else if (it.first=="n")                               // -------------- neutral hadron resolution
                energyResolution=sqrt(pow(0.38/sqrt(jet.E()),2)+pow(0.01,2));
            else if (it.first=="l")                               // ---------------------- lepton resolution
                energyResolution=sqrt(pow(0.00001*jet.pt(),2)+pow(0.005,2));
            else {                                                // ------------------------ everything else
                std::cout << "WARNING: Input jet has unlisted PDG ID!" << std:: endl;
                energyResolution=0; // assume perfect resolution for everything else (?)
            }
        }
    }

    smearDist->SetParameter(0, 1);  // doesn't matter; GetRandom normalizes TF1 
    smearDist->SetParameter(1, 0);
    smearDist->SetParameter(2, energyResolution);
   
    float smearedPt = smearDist->GetRandom();

    jet*=(smearedPt+jet.pt())/jet.pt();
}

// ----------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> discretizeJet(fastjet::PseudoJet jet, 
                                              bool clusterJet) 
{

    static fastjet::JetDefinition   jetDef(fastjet::antikt_algorithm, RPARAM);    
    static fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    static fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );
    static Double_t etaMin = -10,
                    etaMax = 10,
                    etaRes = 0.05,
                    phiMin = 0, 
                    phiMax = 2*TMath::Pi(),
                    phiRes = 0.05; 
    static Int_t etaNBins = TMath::Floor((etaMax - etaMin)/etaRes), 
                 phiNBins = TMath::Floor((phiMax - phiMin)/phiRes);
    
    std::map<TString,std::pair<TH2D*,TH2D*>> resolutionHistos;
    resolutionHistos["n"] = std::pair<TH2D*,TH2D*>(new TH2D("","",etaNBins,etaMin,etaMax,
                                                                  phiNBins,phiMin,phiMax),
                                                   new TH2D("","",etaNBins,etaMin,etaMax,
                                                                  phiNBins,phiMin,phiMax));

    std::vector< fastjet::PseudoJet > newPar;
    
    // collect constituents
    for(int iConst=0; iConst < jet.constituents().size(); iConst++) {
        fastjet::PseudoJet constituent = jet.constituents().at(iConst);
        smearJetPt(constituent);

        for(auto& it: ids) {
            if (resolutionHistos[it.first].first == 0) continue;
            std::vector<int> tidSet = it.second;

            if (std::find(tidSet.begin(), tidSet.end(), constituent.user_index()) != tidSet.end()) {
               resolutionHistos[it.first].first->Fill(constituent.eta(),
                                                      constituent.phi(),
                                                      constituent.E());
               resolutionHistos[it.first].second->Fill(constituent.eta(),
                                                       constituent.phi(),
                                                       constituent.pt());
            }
        }
    }

    // build new constituents
    for(auto& it: resolutionHistos) {
        TH2D* EHisto =it.second.first;
        TH2D* pTHisto=it.second.second;

        // loop over bins
        for(int iEta=1; iEta <= etaNBins; iEta++) {
        for(int iPhi=1; iPhi <= phiNBins; iPhi++) {
            Int_t binNum = EHisto->GetBin(iEta,iPhi);
            Double_t binEta= EHisto->GetBinCenter(binNum);
            Double_t binPhi= EHisto->GetBinCenter(binNum);
            Double_t binE  = EHisto->GetBinContent(binNum);
            Double_t binpT =pTHisto->GetBinContent(binNum);

            Double_t px = binpT*TMath::Cos(binPhi),
                     py = binpT*TMath::Sin(binPhi),
                     pz = binpT*sinh(binEta);

            // reconstruct PseudoJet using pE 4-vector
            // (this is the only way in fastjet unless we assume m=0)
            fastjet::PseudoJet newConst(px,py,pz,binE);
            newConst.set_user_index(ids[it.first][0]); // just steal some id for now
            newPar.push_back(newConst);
        }}

        // prevent memory leaks
        delete it.second.first;
        delete it.second.second;
    }

    // build jet or just return constituents
    if(clusterJet) {
        fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newPar, jetDef, fjAreaDefinition);
        std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));
        return out_jets;
    } else {
        return newPar;
    }
}
