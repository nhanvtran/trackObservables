#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
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

#include "CmdLine.hh"
#include "EventMixer.hh"
#include "PU14.hh"
#include "puppiContainer.hh"

#include "LHEF.h"

float deltaR (float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = std::abs(phi1-phi2); if (dphi>M_PI) dphi-=2*M_PI;  
    return sqrt(deta*deta + dphi*dphi);
}

double findOption(const TString& option, const TString& tag) {
   TObjArray *toks = (TObjArray*) tag.Tokenize(":"); 
   Double_t output=-1;

   for(int i=0; i < toks->GetEntries(); i++) {
        TString tok = ((TObjString*) toks->At(i))->String();

        if(!tok.Contains(option)) continue;
        output=tok.ReplaceAll(option,"").Atof();
        break; 
   }

   std::cout<<"Setting option " << option << " to " << output << std::endl;
   delete toks;
   return output;
}


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
std::vector<float> j_d2_a1_b1;
std::vector<float> j_d2_a1_b2;
std::vector<float> j_n2_b1;
std::vector<float> j_n2_b2;
std::vector<float> j_m2_b1;
std::vector<float> j_m2_b2;

std::vector<float> j_tau1_b1_mmdt;
std::vector<float> j_tau2_b1_mmdt;
std::vector<float> j_tau3_b1_mmdt;
std::vector<float> j_tau1_b2_mmdt;
std::vector<float> j_tau2_b2_mmdt;
std::vector<float> j_tau3_b2_mmdt;
std::vector<float> j_tau21_b2_mmdt;
std::vector<float> j_tau21_b1_mmdt;
std::vector<float> j_tau32_b2_mmdt;
std::vector<float> j_tau32_b1_mmdt;
std::vector<float> j_c1_b0_mmdt;
std::vector<float> j_c1_b1_mmdt;
std::vector<float> j_c1_b2_mmdt;
std::vector<float> j_c2_b1_mmdt;
std::vector<float> j_c2_b2_mmdt;
std::vector<float> j_d2_b1_mmdt;
std::vector<float> j_d2_b2_mmdt;
std::vector<float> j_d2_a1_b1_mmdt;
std::vector<float> j_d2_a1_b2_mmdt;
std::vector<float> j_n2_b1_mmdt;
std::vector<float> j_n2_b2_mmdt;
std::vector<float> j_m2_b1_mmdt;
std::vector<float> j_m2_b2_mmdt;

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

// smearing globals and functions
float resFudgeFactor = 1.;

TRandom3 *smearDist;
void smearJetPt(fastjet::PseudoJet &jet);
std::vector<fastjet::PseudoJet> discretizeEvent(std::vector<fastjet::PseudoJet> &particles, 
                                                Double_t discretize=0, 
                                                Double_t discretizeEcal=0, 
                                                Double_t numberOfPileup=0, 
                                                Double_t maxChargedPt=1e10, 
                                                Double_t maxChargedDr=1e10, 
                                                Double_t trackingEfficiency=1.0);


/////////////////////////////////////////////////////////////////////////////////////////
// MAIN LOOP                                                                           //
/////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char **argv) {
    
    std::string type  = argv[1];   // type "gg" or "qq"
    std::string indir = argv[2];   // where to find input files 
    int min = atoi(argv[3]);      // events to run over
    int max = atoi(argv[4]);      // events to run over
    std::string tag   = argv[5];      // detector type
    std::string jobid = argv[6];      // condorization job id 
    
    RPARAM = 0.4;
    smearDist = new TRandom3();

    char inName[192];
    sprintf( inName, "%s/%s.lhe",indir.c_str(),type.c_str() );
    // sprintf( inName, "/uscms_data/d3/ecoleman/TrackObservablesStudy/trackObservables/processing/pythia82-lhc13-WW-pt1-50k-2.lhe" );
    std::cout << "fname = " << inName << std::endl;
    std::ifstream ifsbkg (inName) ;
    LHEF::Reader reader(ifsbkg) ;

    char outName[192];
    sprintf(outName, "processed-%s-%s.root", type.c_str(), jobid.c_str());
    TFile *f = TFile::Open(outName,"RECREATE");
    f->cd();
    TTree *t_tracks = new TTree("t_tracks","Tree with vectors");
    TTree *t_tragam = new TTree("t_tragam","Tree with vectors");
    TTree *t_allpar = new TTree("t_allpar","Tree with vectors");
    declareBranches(t_tracks);
    declareBranches(t_tragam);
    declareBranches(t_allpar);
  
    std::vector<string> pufiles;
    pufiles.push_back("program");
    pufiles.push_back("-hard");
    pufiles.push_back("lhc14-pythia8-4C-minbias-nev100.pu14.gz");
    pufiles.push_back("-pileup");
    pufiles.push_back("lhc14-pythia8-4C-minbias-nev100.pu14.gz");
    pufiles.push_back("-npu");
    pufiles.push_back("20");
    CmdLine cmdline(pufiles);
    EventMixer* mixer;
    
    if (tag.find("p")!=std::string::npos) {
      mixer=new EventMixer(&cmdline);
    }
          
    // Modify resolution fudge factor if desired
    resFudgeFactor = findOption("r",tag);
    if(resFudgeFactor<=0) resFudgeFactor = 1;

    evtCtr = 0;
    std::vector < fastjet::PseudoJet > particles;
    std::vector < fastjet::PseudoJet > newparticles;
    // loop over events
    while ( reader.readEvent () ) {
        
        ++evtCtr;
        if (evtCtr < min) continue;
        if (evtCtr > max) break;
        
        if (true /*evtCtr % 100 == 0*/) std::cout << "event " << evtCtr << "\n";
        
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
                curpar.set_user_info(new PU14(pdgid,-1,-1));
                particles.push_back( curpar );
            }   

        }

        if (tag.find("p")!=std::string::npos) {

          if (!mixer->next_event()) { // when running out of PU events start from the beginning
            delete mixer;
            mixer=new EventMixer(&cmdline);
            mixer->next_event();
          }

          vector<PseudoJet> full_event = mixer->particles() ;
          vector<PseudoJet> hard_event, pileup_event;
          SelectorIsHard().sift(full_event, hard_event, pileup_event); // this sifts the full event into two vectors

          if (tag.find("i")!=std::string::npos) {
             puppiContainer curEvent(particles, pileup_event);
             particles = curEvent.puppiFetch(20);
             for (unsigned int i = 0; i < particles.size(); ++i)
               particles[i].set_user_index( particles[i].user_info<PU14>().pdg_id() );
          } else {
             for (unsigned int i = 0; i < pileup_event.size(); ++i) {
               if(pileup_event[i].user_info<PU14>().charge()==0) continue;
               pileup_event[i].set_user_index( pileup_event[i].user_info<PU14>().pdg_id() );
               particles.push_back(pileup_event[i]);
             }
          }
        }

        // Discretize neutral hadrons
        static Double_t discretize = findOption("h",tag); //(tag.find("h")!=std::string::npos);
        static Double_t numberOfPileup=(tag.find("q")!=std::string::npos) ? findOption("q",tag) : 0;
        static Double_t discretizeEcal = findOption("e",tag); //(tag.find("e")!=std::string::npos);

        // Threshold above which track reconstruction in jet core is expected 
        // to fail and charged hadrons are reconstructed as neutral hadrons. 
        // Take value where CMS HCAL resolution gets better than tracking resolution
        
        //static Double_t maxChargedPt=(tag.find("t")!=std::string::npos) ? 110 : 1e10;
        static Double_t maxChargedPt=(tag.find("t")!=std::string::npos) ? findOption("t", tag) : 1e10;

        // Distance to nearest neighbor below which track reconstruction in jet 
        // core is expected to fail and charged hadrons are reconstructed as neutral hadrons.
        static Double_t maxChargedDr=(tag.find("s")!=std::string::npos) ? 0.01 : 1e10; 

        // Tracking efficiency
        static Double_t trackingEfficiency=(tag.find("u")!=std::string::npos) ? 0.9 : 1.0; 

        // Discretize if desired
        if (discretize > 0) newparticles = discretizeEvent(particles,
                                                           discretize,
                                                           discretizeEcal,
                                                           numberOfPileup, 
                                                           maxChargedPt, 
                                                           maxChargedDr, 
                                                           trackingEfficiency);
        else newparticles = particles;

        // Smear particle momenta
        if (tag.find("r")!=std::string::npos) 
          for (unsigned int i = 0; i < newparticles.size(); ++i)
            smearJetPt(newparticles[i]);

        //std::cout << "number of particles = " << newparticles.size() 
        //          << ", " << particles.size() 
        //          << ", " << float(newparticles.size())/float(particles.size()) << std::endl;
        analyzeEvent( newparticles, t_tracks, t_tragam, t_allpar );        
    }
    
    std::cout << "finish loop" << std::endl;
    
    f->cd();
    t_tracks->Write();
    t_tragam->Write();
    t_allpar->Write();    
    f->Close();
    
    return 0 ;
}



/////////////////////////////////////////////////////////////////////////////////////////
// MAIN LOOP                                                                           //
/////////////////////////////////////////////////////////////////////////////////////////

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

    if(out_jets.size()>0) {
       thisClustering->delete_self_when_unused();
    } else delete thisClustering;
    
    // std::cout << "delete clustering" << std::endl;

}


/////////////////////////////////////////////////////////////////////////////////////////
// declareBranches:                                                                    // 
// Add branches into the analysis                                                      //
/////////////////////////////////////////////////////////////////////////////////////////

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
    t->Branch("j_d2_a1_b1"       , &j_d2_a1_b1       );
    t->Branch("j_d2_a1_b2"       , &j_d2_a1_b2       );
    t->Branch("j_m2_b1"          , &j_m2_b1          );
    t->Branch("j_m2_b2"          , &j_m2_b2          );
    t->Branch("j_n2_b1"          , &j_n2_b1          );
    t->Branch("j_n2_b2"          , &j_n2_b2          );

    t->Branch("j_tau1_b1_mmdt"        , &j_tau1_b1_mmdt        );
    t->Branch("j_tau2_b1_mmdt"        , &j_tau2_b1_mmdt        );
    t->Branch("j_tau3_b1_mmdt"        , &j_tau3_b1_mmdt        );    
    t->Branch("j_tau1_b2_mmdt"        , &j_tau1_b2_mmdt        );
    t->Branch("j_tau2_b2_mmdt"        , &j_tau2_b2_mmdt        );
    t->Branch("j_tau3_b2_mmdt"        , &j_tau3_b2_mmdt        );    
    t->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    t->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    t->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    t->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    t->Branch("j_tau32_b1_mmdt"       , &j_tau32_b1_mmdt       );
    t->Branch("j_tau32_b2_mmdt"       , &j_tau32_b2_mmdt       );
    t->Branch("j_c1_b0_mmdt"          , &j_c1_b0_mmdt          );
    t->Branch("j_c1_b1_mmdt"          , &j_c1_b1_mmdt          );
    t->Branch("j_c1_b2_mmdt"          , &j_c1_b2_mmdt          );
    t->Branch("j_c2_b1_mmdt"          , &j_c2_b1_mmdt          );
    t->Branch("j_c2_b2_mmdt"          , &j_c2_b2_mmdt          );
    t->Branch("j_d2_b1_mmdt"          , &j_d2_b1_mmdt          );
    t->Branch("j_d2_b2_mmdt"          , &j_d2_b2_mmdt          );
    t->Branch("j_d2_a1_b1_mmdt"       , &j_d2_a1_b1_mmdt       );
    t->Branch("j_d2_a1_b2_mmdt"       , &j_d2_a1_b2_mmdt       );
    t->Branch("j_m2_b1_mmdt"          , &j_m2_b1_mmdt          );
    t->Branch("j_m2_b2_mmdt"          , &j_m2_b2_mmdt          );
    t->Branch("j_n2_b1_mmdt"          , &j_n2_b1_mmdt          );
    t->Branch("j_n2_b2_mmdt"          , &j_n2_b2_mmdt          );

    t->Branch("j_mass_trim"      , &j_mass_trim      );
    t->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    t->Branch("j_mass_prun"      , &j_mass_prun      );
    t->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    t->Branch("j_mass_sdm1"      , &j_mass_sdm1      );
    t->Branch("j_multiplicity"   , &j_multiplicity   );

}

/////////////////////////////////////////////////////////////////////////////////////////
// clearVectors:                                                                       // 
// Cleanup after each tree fill()                                                      //
/////////////////////////////////////////////////////////////////////////////////////////

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
    j_d2_a1_b1.clear();
    j_d2_a1_b2.clear();
    j_m2_b1.clear();
    j_m2_b2.clear();
    j_n2_b1.clear();
    j_n2_b2.clear();

    j_tau1_b1_mmdt.clear();
    j_tau2_b1_mmdt.clear();
    j_tau3_b1_mmdt.clear();    
    j_tau1_b2_mmdt.clear();
    j_tau2_b2_mmdt.clear();
    j_tau3_b2_mmdt.clear();    
    j_tau21_b1_mmdt.clear();
    j_tau21_b2_mmdt.clear();
    j_tau32_b1_mmdt.clear();
    j_tau32_b2_mmdt.clear();
    j_c1_b0_mmdt.clear();
    j_c1_b1_mmdt.clear();
    j_c1_b2_mmdt.clear();
    j_c2_b1_mmdt.clear();
    j_c2_b2_mmdt.clear();
    j_d2_b1_mmdt.clear();
    j_d2_b2_mmdt.clear();
    j_d2_a1_b1_mmdt.clear();
    j_d2_a1_b2_mmdt.clear();
    j_m2_b1_mmdt.clear();
    j_m2_b2_mmdt.clear();
    j_n2_b1_mmdt.clear();
    j_n2_b2_mmdt.clear();

    j_qjetVol.clear();
    j_mass_trim.clear();
    j_mass_mmdt.clear();
    j_mass_prun.clear();
    j_mass_sdb2.clear();
    j_mass_sdm1.clear();
    j_multiplicity.clear();
}


/////////////////////////////////////////////////////////////////////////////////////////
// PushBackJetInformation:                                                             // 
// Run variable computations on a given jet, fill branch arrays                        //
/////////////////////////////////////////////////////////////////////////////////////////

void PushBackJetInformation(fastjet::PseudoJet jet, int particleContentFlag){

    // for reclustering on the fly....
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, RPARAM);    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );

    // std::cout << "old jet pt = " << jet.pt() << ", and particle content flag = " << particleContentFlag << std::endl;

    fastjet::PseudoJet curjet;  
    std::vector< fastjet::PseudoJet > newparticles;
    bool jet_has_no_particles = false;
    if (particleContentFlag == 0) { 
        curjet = jet;
        jet_has_no_particles = curjet[0]==0 && curjet[1]==0 && 
                               curjet[2]==0 && curjet[3]==0;
    }
    else if (particleContentFlag > 0) {

        std::vector< fastjet::PseudoJet > discretePar;
        discretePar = jet.constituents(); 
        
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
        std::vector<fastjet::PseudoJet> out_jets;
        if (newparticles.size() > 0){
            // cluster them
            fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newparticles, jetDef, fjAreaDefinition);
            out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));        
            
            // std::cout << "now fill " << out_jets.size() <<  std::endl;
            // fill into curjet
            if(out_jets.size()>0) {
              curjet = out_jets[0];
              thisClustering->delete_self_when_unused();
            } else delete thisClustering;
        }
        if (out_jets.size()==0) {
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
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b1(2,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b2(2,2,fastjet::contrib::EnergyCorrelator::pt_R);

    fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b1 (1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b2 (2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b1 (1.0,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b2 (1.0,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);

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
        j_c1_b0.push_back( ECF_C1_b0(curjet) );
        j_c1_b1.push_back( ECF_C1_b1(curjet) );
        j_c1_b2.push_back( ECF_C1_b2(curjet) );
        j_c2_b1.push_back( ECF_C2_b1(curjet) );
        j_c2_b2.push_back( ECF_C2_b2(curjet) );
        j_d2_b1.push_back( ECF_D2_b1(curjet) );
        j_d2_b2.push_back( ECF_D2_b2(curjet) );
        j_d2_a1_b1.push_back( ECF_D2_a1_b1(curjet) );
        j_d2_a1_b2.push_back( ECF_D2_a1_b2(curjet) );
        j_m2_b1.push_back( ECF_M2_b1(curjet) );
        j_m2_b2.push_back( ECF_M2_b2(curjet) );
        j_n2_b1.push_back( ECF_N2_b1(curjet) );
        j_n2_b2.push_back( ECF_N2_b2(curjet) );

        // Groomed variables
        fastjet::PseudoJet mmdtjet=soft_drop_mmdt(curjet);

        j_tau1_b1_mmdt.push_back(  nSub1KT_b1(mmdtjet) );        
        j_tau2_b1_mmdt.push_back(  nSub2KT_b1(mmdtjet) );        
        j_tau3_b1_mmdt.push_back(  nSub3KT_b1(mmdtjet) );        
        j_tau1_b2_mmdt.push_back(  nSub1KT_b2(mmdtjet) );        
        j_tau2_b2_mmdt.push_back(  nSub2KT_b2(mmdtjet) );  
        j_tau3_b2_mmdt.push_back(  nSub3KT_b2(mmdtjet) );  
        j_tau21_b1_mmdt.push_back( nSub2KT_b1(mmdtjet) / nSub1KT_b1(mmdtjet) );
        j_tau21_b2_mmdt.push_back( nSub2KT_b2(mmdtjet) / nSub1KT_b2(mmdtjet) );
        j_tau32_b1_mmdt.push_back( nSub3KT_b1(mmdtjet) / nSub2KT_b1(mmdtjet) );
        j_tau32_b2_mmdt.push_back( nSub3KT_b2(mmdtjet) / nSub2KT_b2(mmdtjet) );
        j_c1_b0_mmdt.push_back( ECF_C1_b0(mmdtjet) );
        j_c1_b1_mmdt.push_back( ECF_C1_b1(mmdtjet) );
        j_c1_b2_mmdt.push_back( ECF_C1_b2(mmdtjet) );
        j_c2_b1_mmdt.push_back( ECF_C2_b1(mmdtjet) );
        j_c2_b2_mmdt.push_back( ECF_C2_b2(mmdtjet) );
        j_d2_b1_mmdt.push_back( ECF_D2_b1(mmdtjet) );
        j_d2_b2_mmdt.push_back( ECF_D2_b2(mmdtjet) );
        j_d2_a1_b1_mmdt.push_back( ECF_D2_a1_b1(mmdtjet) );
        j_d2_a1_b2_mmdt.push_back( ECF_D2_a1_b2(mmdtjet) );
        j_m2_b1_mmdt.push_back( ECF_M2_b1(mmdtjet) );
        j_m2_b2_mmdt.push_back( ECF_M2_b2(mmdtjet) );
        j_n2_b1_mmdt.push_back( ECF_N2_b1(mmdtjet) );
        j_n2_b2_mmdt.push_back( ECF_N2_b2(mmdtjet) );
        
        j_mass_trim.push_back( trimmer1( curjet ).m() );
        j_mass_prun.push_back( pruner1( curjet ).m() );    
        j_mass_mmdt.push_back( mmdtjet.m() );
        j_mass_sdb2.push_back( soft_drop_sdb2( curjet ).m() );
        j_mass_sdm1.push_back( soft_drop_sdm1( curjet ).m() );
        
        j_multiplicity.push_back( (float) curjet.constituents().size() );

    }
    else {

        // N-subjettiness
        j_tau1_b1.push_back( -99 );        
        j_tau2_b1.push_back( -99 );        
        j_tau3_b1.push_back( -99 );        
        j_tau1_b2.push_back( -99 );        
        j_tau2_b2.push_back( -99 );  
        j_tau3_b2.push_back( -99 );  
        j_tau21_b1.push_back( -99 );
        j_tau21_b2.push_back( -99 );
        j_tau32_b1.push_back( -99 );
        j_tau32_b2.push_back( -99 );
        
        // energy correlator     
        j_zlogz.push_back( -99 );   
        j_c1_b0.push_back( -99 );
        j_c1_b1.push_back( -99 );
        j_c1_b2.push_back( -99 );
        j_c2_b1.push_back( -99 );
        j_c2_b2.push_back( -99 );
        j_d2_b1.push_back( -99 );
        j_d2_b2.push_back( -99 );
        j_d2_a1_b1.push_back( -99 );
        j_d2_a1_b2.push_back( -99 );
        j_m2_b1.push_back( -99 );
        j_m2_b2.push_back( -99 );
        j_n2_b1.push_back( -99 );
        j_n2_b2.push_back( -99 );

        // Groomed variables
        j_tau1_b1_mmdt.push_back( -99 );        
        j_tau2_b1_mmdt.push_back( -99 );        
        j_tau3_b1_mmdt.push_back( -99 );        
        j_tau1_b2_mmdt.push_back( -99 );	     
        j_tau2_b2_mmdt.push_back( -99 );  
        j_tau3_b2_mmdt.push_back( -99 );  
        j_tau21_b1_mmdt.push_back( -99 );
        j_tau21_b2_mmdt.push_back( -99 );
        j_tau32_b1_mmdt.push_back( -99 );
        j_tau32_b2_mmdt.push_back( -99 );
        j_c1_b0_mmdt.push_back( -99 );
        j_c1_b1_mmdt.push_back( -99 );
        j_c1_b2_mmdt.push_back( -99 );
        j_c2_b1_mmdt.push_back( -99 );
        j_c2_b2_mmdt.push_back( -99 );
        j_d2_b1_mmdt.push_back( -99 );
        j_d2_b2_mmdt.push_back( -99 );
        j_d2_a1_b1_mmdt.push_back( -99 );
        j_d2_a1_b2_mmdt.push_back( -99 );
        j_m2_b1_mmdt.push_back( -99 );
        j_m2_b2_mmdt.push_back( -99 );
        j_n2_b1_mmdt.push_back( -99 );
        j_n2_b2_mmdt.push_back( -99 );
        
        j_mass_trim.push_back( -99 );
        j_mass_prun.push_back( -99 );
        j_mass_mmdt.push_back( -99 );
        j_mass_sdb2.push_back( -99 );
        j_mass_sdm1.push_back( -99 );
        
        j_multiplicity.push_back( -99 );

    }
}


/////////////////////////////////////////////////////////////////////////////////////////
// smearJetPt                                                                          // 
// Apply resolution smearing to a given pseudojet                                      //
//  - applied resolution depends on PDG ID of jet                                      //
/////////////////////////////////////////////////////////////////////////////////////////

void smearJetPt(fastjet::PseudoJet &jet) {

    Double_t energyResolution=0;
    for(auto& it: ids) { 
        std::vector<int> tidSet = it.second;

        if (std::find(tidSet.begin(), tidSet.end(), jet.user_index()) != tidSet.end()) {
            if(it.first=="c") {                                    // -------------- charged hadron resolution
                //energyResolution=sqrt(pow(0.0001*jet.pt(),2)+pow(0.005,2)); // CMS TDR

                // Delphes CMS tuning https://github.com/sethzenz/Delphes/blob/master/Cards/CMS_Phase_I_NoPileUp.tcl#L159
                energyResolution=sqrt(pow(0.00025*jet.pt(),2)+pow(0.015,2)); 
            } else if (it.first=="p") {                              // ---------------------- photon resolution
                // CMS TDR
                //energyResolution=sqrt(pow(0.027/sqrt(jet.e()),2)+pow(0.15/jet.e(),2)+pow(0.005,2)); 

                // CMS TDR + factor 1.6=1+0.6 for the fact that 60% charged hadrons (leaving 30% in ECAL) are substracted before 30% photons are reconstructed
                //energyResolution=sqrt(pow(0.027/sqrt(1.6*jet.e()),2)+pow(0.15/1.6/jet.e(),2)+pow(1.6*0.005,2)); 

                // CMS TDR + factor 1.6=1+0.6 for the fact that 60% charged hadrons (leaving 30% in ECAL) are substracted before 30% photons are reconstructed
                energyResolution=sqrt(pow(0.027/sqrt(jet.e()/1.6),2)+pow(0.15/jet.e()*1.6,2)+pow(0.005,2)); 

                // CMS tuning https://github.com/cms-met/cmssw/blob/b742cc16aff1915b275cf0847dcff93aa6deab14/RecoMET/METProducers/python/METSigParams_cfi.py#L37
                //energyResolution=sqrt(pow(0.042/jet.e(),2)+pow(0.1/sqrt(jet.e()),2)+pow(0.005,2)); 
            } else if (it.first=="n") {                              // -------------- neutral hadron resolution

                // CMS JET JINST 
                //energyResolution=sqrt(pow(1.20/sqrt(jet.e()),2)+pow(0.05,2)); 

                // CMS JET JINST + factor 7=1+6 for the fact that 60% charged hadrons are substracted before 10% neutral hadrons are reconstructed
                //energyResolution=sqrt(pow(1.20/sqrt(7.*jet.e()),2)+pow(7.*0.05,2)); 

                // CMS JET JINST + factor 7=1+6 for the fact that 60% charged hadrons are substracted before 10% neutral hadrons are reconstructed
                energyResolution=sqrt(pow(1.20/sqrt(jet.e()/7.),2)+pow(0.05,2)); 

                // CMS tuning https://github.com/cms-met/cmssw/blob/b742cc16aff1915b275cf0847dcff93aa6deab14/RecoMET/METProducers/python/METSigParams_cfi.py#L37
                //energyResolution=sqrt(pow(0.41/jet.e(),2)+pow(0.52/sqrt(jet.e()),2)+pow(0.25,2)); 

            } else if (it.first=="l") {                              // ---------------------- lepton resolution
                // CMS TDR
                //energyResolution=sqrt(pow(0.0001*jet.pt(),2)+pow(0.005,2)); 

                // Delphes CMS tuning https://github.com/sethzenz/Delphes/blob/master/Cards/CMS_Phase_I_NoPileUp.tcl#L159
                energyResolution=sqrt(pow(0.00025*jet.pt(),2)+pow(0.015,2)); 
            } else {                                                // ------------------------ everything else
                std::cout << "WARNING: Input jet has unlisted PDG ID!" << std:: endl;
                energyResolution=0; // assume perfect resolution for everything else (?)
            }

            // std::cout << it.first << " " << energyResolution << std::endl;
        }
    }

    Double_t smearedPt = std::max(1e-10,smearDist->Gaus(1,energyResolution*resFudgeFactor));

    jet.reset_momentum(jet.px() * smearedPt,
                       jet.py() * smearedPt,
                       jet.pz() * smearedPt,
                       jet.e()  * smearedPt);

}

/////////////////////////////////////////////////////////////////////////////////////////
// discretizeEvent                                                                     // 
// Apply calorimeter discretization to an event                                        //
//  - option to toggle ecal discretization separately                                  //
/////////////////////////////////////////////////////////////////////////////////////////

std::vector<fastjet::PseudoJet> discretizeEvent(std::vector<fastjet::PseudoJet> &particles, 
                                                Double_t discretize,
                                                Double_t discretizeEcal, 
                                                Double_t numberOfPileup, 
                                                Double_t maxChargedPt, 
                                                Double_t maxChargedDr, 
                                                Double_t trackingEfficiency)
{
    if(discretize<=0) return particles;

    static Double_t etaMin = -3,
                    etaMax = 3,
                    //etaRes = 0.087/sqrt(12), // CMS TDR, factor 3 since particle flow cluster algorithm reaches this precision
                    etaRes = discretize, 
                    phiMin = 0, 
                    phiMax = 2*TMath::Pi(),
                    //phiRes = 0.087/sqrt(12).; // CMS TDR, factor 3 since particle flow cluster algorithm reaches this precision
                    phiRes = discretize; 
    static Int_t etaNBins = TMath::Floor((etaMax - etaMin)/etaRes), 
                 phiNBins = TMath::Floor((phiMax - phiMin)/phiRes);

    TH2D* hcalGrid = new TH2D("hcalGrid","hcalGrid",etaNBins,etaMin,etaMax,phiNBins,phiMin,phiMax);                   

    //static Double_t etaResEcal = 0.017, // CMS ECAL JINST
    //                phiResEcal = 0.017; // CMS ECAL JINST
    static Double_t etaResEcal = (discretizeEcal > 0 ? discretizeEcal : 1),
                    phiResEcal = (discretizeEcal > 0 ? discretizeEcal : 1);
    static Int_t etaNBinsEcal = TMath::Floor((etaMax - etaMin)/etaResEcal), // CMS ECAL JINST
                 phiNBinsEcal = TMath::Floor((phiMax - phiMin)/phiResEcal); // CMS ECAL JINST

    TH2D* ecalGrid = new TH2D("ecalGrid","ecalGrid",etaNBinsEcal,etaMin,etaMax,phiNBinsEcal,phiMin,phiMax);  
    
    // Add in neutral pileup: 0.3 * "q" GeV total in detector on average per event
    static Double_t pileupEnergyInDetector = numberOfPileup * 0.3 ; 
    double addedPileupEnergy = 0;
    while(addedPileupEnergy < pileupEnergyInDetector) {
        int randPhi = smearDist->Uniform(0,2*TMath::Pi());
        int randEta = smearDist->Uniform(-3,3);
        double randEnergy = smearDist->Exp(1/6)+0.5; // e^-6t only above 500 MeV

        ecalGrid->Fill(randEta,randPhi,randEnergy);
        hcalGrid->Fill(randEta,randPhi,randEnergy);
    
        addedPileupEnergy += randEnergy;
    }

    std::vector<fastjet::PseudoJet> newparticles;
    for (unsigned int i = 0; i < particles.size(); ++i){
        // newparticles.push_back( fastjet::PseudoJet(particles[i]) );
        for(auto& it: ids) { 
            std::vector<int> tidSet = it.second;

            if (std::find(tidSet.begin(), tidSet.end(), particles[i].user_index()) != tidSet.end()) {
                bool track=(it.first=="c");
                if(track&&(trackingEfficiency!=1.0)&&(smearDist->Integer(100)>trackingEfficiency*100.0))
                  track=false;
                if(track&&(maxChargedDr<10)) {
                  double cutDr=smearDist->Integer(100)/50.*maxChargedDr;
                  for (unsigned int j = 0; j < particles.size(); ++j) {
                    if((i!=j)&&(deltaR(particles[i].eta(),particles[i].phi(),particles[j].eta(),particles[j].phi())<cutDr)) {
                      track=false;
                      break;
                    }
                  }
                }
                if(track && (maxChargedPt<10000) && (particles[i].pt()>maxChargedPt))
                   track=false;
                if( track || it.first=="l" || (discretizeEcal <= 0 && it.first=="p"))
                {                                   
                    newparticles.push_back( fastjet::PseudoJet(particles[i]) );
                    // std::cout << "c par id = " << particles[i].user_index() << std::endl;
                }
                else if (it.first=="p"){
                    double curphi = particles[i].phi();
                    double cureta = particles[i].eta();
                    int ieta = ecalGrid->GetXaxis()->FindBin( cureta );
                    int iphi = ecalGrid->GetYaxis()->FindBin( curphi );
                    ecalGrid->SetBinContent( ieta,iphi, ecalGrid->GetBinContent(ieta,iphi)+particles[i].e() );
                    // std::cout << "n par id = " << particles[i].user_index() 
                    //                            << "," << particles[i].e() 
                    //                            << "," << ieta << "," 
                    //                            << iphi << "," 
                    //                            << cureta << "," 
                    //                            << curphi 
                    //                            << std::endl;
                }
                else if ((it.first=="n")||(!track)){
                    double curphi = particles[i].phi();
                    double cureta = particles[i].eta();
                    int ieta = hcalGrid->GetXaxis()->FindBin( cureta );
                    int iphi = hcalGrid->GetYaxis()->FindBin( curphi );
                    hcalGrid->SetBinContent( ieta,iphi, hcalGrid->GetBinContent(ieta,iphi)+particles[i].e() );
                    // std::cout << "n par id = " << particles[i].user_index() 
                    //                            << "," << particles[i].e() 
                    //                            << "," << ieta << "," 
                    //                            << iphi << "," 
                    //                            << cureta << "," 
                    //                            << curphi 
                    //                            << std::endl;
                }
                else {                                                // ------------------------ everything else
                    std::cout << "WARNING: Input particle has unlisted PDG ID!" << std:: endl;
                }
            }
        }
    }

    int hcalcellctr = 0;
    for (int i = 0; i < etaNBins; ++i){
        for (int j = 0; j < phiNBins; ++j){
            float curbincontent = hcalGrid->GetBinContent(i+1,j+1);
            if (curbincontent > 0){ 
                float celleta = hcalGrid->GetXaxis()->GetBinCenter(i+1);
                float cellphi = hcalGrid->GetYaxis()->GetBinCenter(j+1);
                float celle   = hcalGrid->GetBinContent(i+1,j+1);
                // float cellpt  = sqrt(celle*celle*(2*exp(2*celleta))/(1+exp(2*celleta)));
                // float cellpt  = sqrt(2*celle*celle/(1+exp(2*celleta)));
                float cellpt = celle*2/(exp(celleta)+exp(-celleta));
                fastjet::PseudoJet curcell = fastjet::PseudoJet(0,0,0,0);
                curcell.reset_PtYPhiM(cellpt,celleta,cellphi,0.0);
                // std::cout << "celle = "    << celle 
                //           << ", cellpt = " << cellpt 
                //           << ", orig e = " << curbincontent 
                //           << ", cell eta = " << celleta << std::endl;
                curcell.set_user_index( 130 );
                newparticles.push_back(curcell);
                hcalcellctr++;
            }
            // std::cout << i+1 << "," << j+1 << "," << curbincontent << std::endl;
        }
    }

    // std::cout << "hcalcellctr = " << hcalcellctr << std::endl;

    int ecalcellctr = 0;
    for (int i = 0; i < etaNBinsEcal; ++i){
        for (int j = 0; j < phiNBinsEcal; ++j){
            float curbincontent = ecalGrid->GetBinContent(i+1,j+1);
            if (curbincontent > 0){ 
                float celleta = ecalGrid->GetXaxis()->GetBinCenter(i+1);
                float cellphi = ecalGrid->GetYaxis()->GetBinCenter(j+1);
                float celle   = ecalGrid->GetBinContent(i+1,j+1);
                // float cellpt  = sqrt(celle*celle*(2*exp(2*celleta))/(1+exp(2*celleta)));
                // float cellpt  = sqrt(2*celle*celle/(1+exp(2*celleta)));
                float cellpt = celle*2/(exp(celleta)+exp(-celleta));
                fastjet::PseudoJet curcell = fastjet::PseudoJet(0,0,0,0);
                curcell.reset_PtYPhiM(cellpt,celleta,cellphi,0.0);
                // std::cout << "celle = "    << celle 
                //           << ", cellpt = " << cellpt 
                //           << ", orig e = " << curbincontent 
                //           << ", cell eta = " << celleta << std::endl;
                curcell.set_user_index( 111 );
                newparticles.push_back(curcell);
                ecalcellctr++;
            }
            // std::cout << i+1 << "," << j+1 << "," << curbincontent << std::endl;
        }
    }

    // std::cout << "ecalcellctr = " << ecalcellctr << std::endl;

    delete hcalGrid;
    delete ecalGrid;
    return newparticles;
}
