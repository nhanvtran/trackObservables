#include "HepMC/IO_GenEvent.h"
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
#include "RecursiveSoftDrop.hh"

#include "EnergyCorrelations.hh"
#include "CmdLine.hh"
#include "EventMixer.hh"
#include "PU14.hh"
#include "puppiContainer.hh"

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

std::string intToString(int i)
{
    std::stringstream ss;
    std::string s;
    ss << i;
    s = ss.str();

    return s;
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
EnergyCorrelations* fECF;
int evtCtr;
float RPARAM;

int njets;
std::vector<int> j_nbHadrons;
std::vector<int> j_ncHadrons;
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
std::vector<float> j_n3_b1;
std::vector<float> j_n3_b2;
std::vector<float> j_m3_b1;
std::vector<float> j_m3_b2;

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
std::vector<float> j_n3_b1_mmdt;
std::vector<float> j_n3_b2_mmdt;
std::vector<float> j_m3_b1_mmdt;
std::vector<float> j_m3_b2_mmdt;

std::vector<float> j_qjetVol;
std::vector<float> j_mass_trim;
std::vector<float> j_mass_mmdt;
std::vector<float> j_mass_rsd;
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
    { "l" , {-11,11,-13,13,-15,15}                                          },
    { "b" , {511,-511,521,-521,531,-531,541,-541,10513,-10513,10523,-10523,10533,-10533,10543,-10543,20513,-20513,20523,-20523,20533,-20533,20543,-20543,5122,-5122,5132,5112,-5112,5114,-5114,-5132,5142,-5142,5212,-5212,5214,-5214,5222,-5222,5224,-5224,5232,-5232,5242,-5242,5312,-5312,5314,-5314,5322,-5322,5324,-5324,5332,-5332,5334,-5334,5342,-5342,5412,-5412,5414,-5414,5422,-5422,5424,-5424,5432,-5432,5434,-5434,5442,-5442,5444,-5444,5512,-5512,5514,-5514,5522,-5522,5524,-5524,5532,-5532,5534,-5534,5542,-5542,5544,-5544,5554,-5554,513,-513}},
    { "chad", {411,-411,10411,-10411,10413,-10413,20413,-20413,421,-421,10421,-10421,10423,-10423,20423,-20423,431,-431,10433,-10433,20433,-20433,4122,-4122,4112,-4112,4212,-4212,4222,-4222,4114,-4114,4214,-4214,4224,-4224,4232,-4232,4132,-4132,4322,-4322,4312,-4312,4324,-4324,4314,-4314,4332,-4332,4334,-4334,4412,-4412,4422,-4422,4414,-4414,4424,-4424,4432,-4432,4434,-4434,4444,-4444}}
};

std::vector<int> c_ids = ids["c"];
std::vector<int> p_ids = ids["p"]; 
std::vector<int> n_ids = ids["n"]; 
std::vector<int> l_ids = ids["l"]; 
std::vector<int> b_ids = ids["b"];
std::vector<int> chad_ids = ids["chad"];

int activeAreaRepeats = 1;
double ghostArea = 0.01;
double ghostEtaMax = 7.0;
    
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, TTree* t_tracks, TTree* t_tragam, TTree* t_allpar, std::vector < fastjet::PseudoJet > HFparticles);
void declareBranches(TTree* t);
void clearVectors();
void PushBackJetInformation(fastjet::PseudoJet jet, int particleContentFlag, std::vector < fastjet::PseudoJet > HFparticles);

// smearing functions
TRandom3 *smearDist;
void smearJetPt(fastjet::PseudoJet &jet, Double_t shiftChargedScale, Double_t shiftPhotonScale, Double_t shiftHadronScale);
std::vector<fastjet::PseudoJet> discretizeJet(fastjet::PseudoJet jet, bool clusterJet = false);
std::vector<fastjet::PseudoJet> discretizeEvent(std::vector<fastjet::PseudoJet> &particles, bool discretizeEcal=false, Double_t numberOfPileup=0, Double_t maxChargedPt=1e10, Double_t maxChargedDr=1e10, Double_t trackingEfficiency=1.0);

////////////////////-----------------------------------------------

int main (int argc, char **argv) {
    
    std::string type  = argv[1];   // input filename
    std::string indir = argv[2];   // where to find input files 
    int min = atoi(argv[3]);      // events to run over
    int max = atoi(argv[4]);      // events to run over
    std::string tag = argv[5];      // detector type
    
    RPARAM = 0.8;
    smearDist = new TRandom3();

    char inName[192];
    bool isLHE = type.find("lhe") != std::string::npos;
    bool isROOT = type.find("root") != std::string::npos;
    sprintf( inName, "%s/%s",indir.c_str(),type.c_str() );
    // sprintf( inName, "/uscms_data/d3/ecoleman/TrackObservablesStudy/trackObservables/processing/pythia82-lhc13-WW-pt1-50k-2.lhe" );
    std::cout << "fname = " << inName << std::endl;
    std::ifstream ifsbkg (inName) ;
    
    char outName[192];
    sprintf(outName, "processed-%s-%s.root", type.c_str(), tag.c_str());
    TFile *f = TFile::Open(outName,"RECREATE");
    f->cd();
    TTree *t_tracks = new TTree("t_tracks","Tree with vectors");
    TTree *t_tragam = new TTree("t_tragam","Tree with vectors");
    TTree *t_allpar = new TTree("t_allpar","Tree with vectors");
    declareBranches(t_tracks);
    declareBranches(t_tragam);
    declareBranches(t_allpar);

    static Double_t numberOfPileup=findOption("p",tag);
  
    std::vector<string> pufiles;
    pufiles.push_back("program");
    pufiles.push_back("-hard");
    pufiles.push_back("lhc14-pythia8-4C-minbias-nev100.pu14.gz");
    pufiles.push_back("-pileup");
    pufiles.push_back("lhc14-pythia8-4C-minbias-nev100.pu14.gz");
    pufiles.push_back("-npu");
    pufiles.push_back(intToString(numberOfPileup));
    CmdLine cmdline(pufiles);
    EventMixer* mixer;
    
    if (tag.find("p")!=std::string::npos) {
      mixer=new EventMixer(&cmdline);
    }
    fECF = new EnergyCorrelations();
    // evtCtr = 0;
    std::vector < fastjet::PseudoJet > particles;
    std::vector < fastjet::PseudoJet > hfparticles;
    std::vector < fastjet::PseudoJet > newparticles;
    std::vector < fastjet::PseudoJet > neutrinos;

    // loop over events
    bool nextEvent = true;

    LHEF::Reader *reader;
    if(isLHE) { 
      reader = new LHEF::Reader(ifsbkg); 
      nextEvent = reader->readEvent();
    }

    HepMC::GenEvent* evt;
    HepMC::IO_GenEvent *hepmcreader;
    if(!isLHE && !isROOT) {
      hepmcreader = new HepMC::IO_GenEvent(inName,std::ios::in);
      evt = hepmcreader->read_next_event();
      nextEvent = (evt != 0);
    }


    TFile *fin;
    TTree *tin;
    const Int_t kMaxEvent = 1;
    const Int_t kMaxGenParticle = 9994;
    Int_t           Event_;
    UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
    UInt_t          Event_fBits[kMaxEvent];   //[Event_]
    Long64_t        Event_Number[kMaxEvent];   //[Event_]
    Int_t           Event_size;    
    Int_t           GenParticle_;
    UInt_t          GenParticle_fUniqueID[kMaxGenParticle];   //[GenParticle_]
    UInt_t          GenParticle_fBits[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_PID[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_Status[kMaxGenParticle];   //[GenParticle_]
    Double_t        GenParticle_E[kMaxGenParticle];   //[GenParticle_]
    Double_t        GenParticle_Px[kMaxGenParticle];   //[GenParticle_]
    Double_t        GenParticle_Py[kMaxGenParticle];   //[GenParticle_]
    Double_t        GenParticle_Pz[kMaxGenParticle];   //[GenParticle_]
    Double_t        GenParticle_PT[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_M1[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_M2[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_D1[kMaxGenParticle];   //[GenParticle_]
    Int_t           GenParticle_D2[kMaxGenParticle];   //[GenParticle_]

    if(isROOT){
        fin = new TFile(inName,"READ");     // Open root input file from madgraph+pythia
        tin = (TTree*)fin->Get("STDHEP");    // Get tree from root file
        tin->SetMakeClass(1); // Necessary here !!!
        tin->SetBranchAddress("Event",                &Event_);
        tin->SetBranchAddress("Event.fUniqueID",       Event_fUniqueID);
        tin->SetBranchAddress("Event.fBits",           Event_fBits);
        tin->SetBranchAddress("Event.Number",          Event_Number);
        tin->SetBranchAddress("Event_size",           &Event_size);
        tin->SetBranchAddress("GenParticle",          &GenParticle_);
        tin->SetBranchAddress("GenParticle.fUniqueID", GenParticle_fUniqueID);
        tin->SetBranchAddress("GenParticle.fBits",     GenParticle_fBits);
        tin->SetBranchAddress("GenParticle.PID",       GenParticle_PID);
        tin->SetBranchAddress("GenParticle.Status",    GenParticle_Status);
        tin->SetBranchAddress("GenParticle.E",         GenParticle_E);
        tin->SetBranchAddress("GenParticle.Px",        GenParticle_Px);
        tin->SetBranchAddress("GenParticle.Py",        GenParticle_Py);
        tin->SetBranchAddress("GenParticle.Pz",        GenParticle_Pz);        
        tin->SetBranchAddress("GenParticle.M1",        GenParticle_M1);
        tin->SetBranchAddress("GenParticle.M2",        GenParticle_M2);
        tin->SetBranchAddress("GenParticle.D1",        GenParticle_D1);
        tin->SetBranchAddress("GenParticle.D2",        GenParticle_D2);
    }

    Long64_t nentries = 0;
    if (isROOT){ 
        nentries = tin->GetEntriesFast();
        nextEvent = true;
    }
    std::cout << "Number of entries = " << nentries << ", " << isROOT << ", " << nextEvent << std::endl;



    // Start event loop
    while ( nextEvent ) {
      
        ++evtCtr;
        if (evtCtr < min) continue;
        if (evtCtr > max) break;      
        if (evtCtr % 100 == 0) std::cout << "event " << evtCtr << "\n";
        
        // per event
        hfparticles.clear();
        particles.clear();
        neutrinos.clear();
        bool onlyWplus = (tag.find("W+")!=std::string::npos);
        bool onlyWminus = (tag.find("W-")!=std::string::npos);
        bool containsWplus=false;
        bool containsWminus=false;
        if(isLHE) { //Read an LHE file
            //std::cout << "reader.hepeup.IDUP.size() = " << reader.hepeup.IDUP.size() << std::endl;
            for (unsigned int i = 0 ; i < reader->hepeup.IDUP.size(); ++i){
                if(reader->hepeup.IDUP.at(i)==24) containsWplus=true;
                if(reader->hepeup.IDUP.at(i)==-24) containsWminus=true;
                //std::cout << reader->hepeup.IDUP.at(i) << " -- " << reader->hepeup.PUP.at(i).at(3) << " -- " << reader->hepeup.ISTUP.at(i) << std::endl;
                if((reader->hepeup.ISTUP.at(i) == 1) && (TMath::Abs(reader->hepeup.IDUP.at(i))!=14) && (TMath::Abs(reader->hepeup.IDUP.at(i))!=16) && (TMath::Abs(reader->hepeup.IDUP.at(i))!=18)){
                    float px = reader->hepeup.PUP.at(i).at(0);
                    float py = reader->hepeup.PUP.at(i).at(1);
                    float pz = reader->hepeup.PUP.at(i).at(2);
                    float e  = reader->hepeup.PUP.at(i).at(3);                                    
                    fastjet::PseudoJet curpar   = fastjet::PseudoJet( px, py, pz, e );
                    int pdgid = reader->hepeup.IDUP.at(i);
                    curpar.set_user_index( pdgid );
                    curpar.set_user_info(new PU14(pdgid,-1,-1));
                    particles.push_back( curpar );
                }                 
            } 
        } 
        else if (isROOT){
            tin->GetEntry(evtCtr);
            //std::cout << "Evt ctr = " << evtCtr << ", particles " << GenParticle_ << std::endl;
            for (int ipar = 0; ipar < GenParticle_; ipar++){
                if(GenParticle_PID[ipar]==24) containsWplus=true;
                if(GenParticle_PID[ipar]==-24) containsWminus=true;
                
                if (GenParticle_Status[ipar] == 1){
                    float px = GenParticle_Px[ipar];
                    float py = GenParticle_Py[ipar];
                    float pz = GenParticle_Pz[ipar];
                    float e  = GenParticle_E[ipar];    
                    fastjet::PseudoJet curpar   = fastjet::PseudoJet( px, py, pz, e );
                    int pdgid = GenParticle_PID[ipar];
                    curpar.set_user_index( pdgid );
                    curpar.set_user_info(new PU14(pdgid,-1,-1));
                    if(TMath::Abs(pdgid) != 14 && TMath::Abs(pdgid) != 16 && TMath::Abs(pdgid) != 18 ){
                        particles.push_back( curpar );                   
                        if (std::find(b_ids.begin(), b_ids.end(), curpar.user_index()) != b_ids.end()) {
                            hfparticles.push_back( curpar );
                        }
                        if (std::find(chad_ids.begin(), chad_ids.end(), curpar.user_index()) != chad_ids.end()) {
                            if(abs(GenParticle_Status[ipar])<=2) hfparticles.push_back( curpar );
                        }                   
                    } else {
                        neutrinos.push_back(curpar);
                    }
                } 
            }
        }
        else { //Read a HepMC file
            for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){
                if((*p)->pdg_id()==24) containsWplus=true;
                if((*p)->pdg_id()==-24) containsWminus=true;
                //std::cout << (*p)->pdg_id() << " -- " << (*p)->momentum().perp() << " -- " << (*p)->status() << std::endl;
                float px = (*p)->momentum().px();
                float py = (*p)->momentum().py();
                float pz = (*p)->momentum().pz();
                float e  = (*p)->momentum().e();
                fastjet::PseudoJet curpar   = fastjet::PseudoJet( px, py, pz, e );
                int pdgid = (*p)->pdg_id();
                if(TMath::Abs(pdgid) != 14 && TMath::Abs(pdgid) != 16 && TMath::Abs(pdgid) != 18 ){
                    curpar.set_user_index( pdgid );
                    curpar.set_user_info(new PU14(pdgid,-1,-1));
                    if (std::find(b_ids.begin(), b_ids.end(), curpar.user_index()) != b_ids.end()) {
                        //if(abs((*p)->status())<=2)
                        hfparticles.push_back( curpar );
                    }
                    if (std::find(chad_ids.begin(), chad_ids.end(), curpar.user_index()) != chad_ids.end()) {
                        if(abs((*p)->status())<=2) hfparticles.push_back( curpar );
                    }
                    if(abs((*p)->status())==1) particles.push_back( curpar );
                }
            }
        }
        
        if (numberOfPileup>0) {
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
                particles = curEvent.puppiFetch(numberOfPileup);
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
      	
        if(((!onlyWplus)||(containsWplus))&&((!onlyWminus)||(containsWminus))) {
            // discretize neutral hadrons
            bool discretize = (tag.find("h")!=std::string::npos);
            static Double_t nPU=numberOfPileup*(tag.find("q")!=std::string::npos);
            static bool discretizeEcal=(tag.find("e")!=std::string::npos);
            static Double_t maxChargedPt=(tag.find("t")!=std::string::npos)?110:1e10; // Threshold above which track reconstruction in jet core is expected to fail and charged hadrons are reconstructed as neutral hadrons. Take value where CMS HCAL resolution gets better than tracking resolution
            static Double_t maxChargedDr=(tag.find("s")!=std::string::npos)?0.01:1e10; // Distance to nearest neighbor below which track reconstruction in jet core is expected to fail and charged hadrons are reconstructed as neutral hadrons.
            static Double_t trackingEfficiency=(tag.find("u")!=std::string::npos)?0.9:1.0; // Tracking efficiency
            if (discretize) newparticles = discretizeEvent(particles,discretizeEcal,nPU, maxChargedPt, maxChargedDr, trackingEfficiency);
            else newparticles = particles;
            // smear particle momenta
            static Double_t shiftChargedScale=findOption("CS",tag);
            static Double_t shiftPhotonScale=findOption("PS",tag);
            static Double_t shiftHadronScale=findOption("HS",tag);
            if (tag.find("r")!=std::string::npos)
                for (unsigned int i = 0; i < newparticles.size(); ++i)
                    smearJetPt(newparticles[i],shiftChargedScale,shiftPhotonScale,shiftHadronScale);
            // std::cout << "number of particles = " << newparticles.size() << ", " << particles.size() << ", " << float(newparticles.size())/float(particles.size()) << std::endl;
            analyzeEvent( newparticles, t_tracks, t_tragam, t_allpar, hfparticles );        
        }
        
        if(isLHE) { 
            nextEvent = reader->readEvent();
        }
        else if(isROOT){
            if (evtCtr > nentries) nextEvent = false;
        }
        else { 
            delete evt;
            *hepmcreader >> evt;
            nextEvent = (evt != 0);
        }
    }
    
    std::cout << "finish loop" << std::endl;
    
    f->cd();
    t_tracks->Write();
    t_tragam->Write();
    t_allpar->Write();    
    f->Close();
    
    return 0 ;
}

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, TTree* t_tracks, TTree* t_tragam, TTree* t_allpar, std::vector < fastjet::PseudoJet > HFparticles){
    
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
        PushBackJetInformation(out_jets[i],0,HFparticles);
    }
    t_allpar->Fill();
    clearVectors();

    for (unsigned int i = 0; i < out_jets.size(); i++){
        PushBackJetInformation(out_jets[i],1,HFparticles);
    }
    t_tracks->Fill();
    clearVectors();

    // std::cout << "tracks" << std::endl;

    for (unsigned int i = 0; i < out_jets.size(); i++){
        PushBackJetInformation(out_jets[i],2,HFparticles);
    }
    t_tragam->Fill();
    clearVectors();
    
    // std::cout << "tracks+photons" << std::endl;

    if(out_jets.size()>0) {
       thisClustering->delete_self_when_unused();
    } else delete thisClustering;
    
    // std::cout << "delete clustering" << std::endl;

}

// ----------------------------------------------------------------------------------
void declareBranches( TTree* t ){

    t->Branch("njets"            , &njets            );
    t->Branch("j_ptfrac"         , &j_ptfrac         );
    t->Branch("j_nbHadrons"	 , &j_nbHadrons      );
    t->Branch("j_ncHadrons"      , &j_ncHadrons      );
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
    t->Branch("j_m3_b1"          , &j_m3_b1          );
    t->Branch("j_m3_b2"          , &j_m3_b2          );
    t->Branch("j_n3_b1"          , &j_n3_b1          );
    t->Branch("j_n3_b2"          , &j_n3_b2          );

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
    t->Branch("j_m3_b1_mmdt"          , &j_m3_b1_mmdt          );
    t->Branch("j_m3_b2_mmdt"          , &j_m3_b2_mmdt          );
    t->Branch("j_n3_b1_mmdt"          , &j_n3_b1_mmdt          );
    t->Branch("j_n3_b2_mmdt"          , &j_n3_b2_mmdt          );

    t->Branch("j_mass_trim"      , &j_mass_trim      );
    t->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    t->Branch("j_mass_rsd"       , &j_mass_rsd      );
    t->Branch("j_mass_prun"      , &j_mass_prun      );
    t->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    t->Branch("j_mass_sdm1"      , &j_mass_sdm1      );
    t->Branch("j_multiplicity"   , &j_multiplicity   );

}

// ----------------------------------------------------------------------------------
void clearVectors(){
    j_pt.clear();
    j_ptfrac.clear();
    j_nbHadrons.clear();
    j_ncHadrons.clear();
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
    j_m3_b1.clear();
    j_m3_b2.clear();
    j_n3_b1.clear();
    j_n3_b2.clear();

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
    j_m3_b1_mmdt.clear();
    j_m3_b2_mmdt.clear();
    j_n3_b1_mmdt.clear();
    j_n3_b2_mmdt.clear();

    j_qjetVol.clear();
    j_mass_trim.clear();
    j_mass_mmdt.clear();
    j_mass_rsd.clear();
    j_mass_prun.clear();
    j_mass_sdb2.clear();
    j_mass_sdm1.clear();
    j_multiplicity.clear();
}

// ----------------------------------------------------------------------------------
void PushBackJetInformation(fastjet::PseudoJet jet, int particleContentFlag, std::vector < fastjet::PseudoJet > HFparticles){

    // for reclustering on the fly....
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, RPARAM);    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );

    // std::cout << "old jet pt = " << jet.pt() << ", and particle content flag = " << particleContentFlag << std::endl;

    fastjet::PseudoJet curjet;  
    int nbHadrons =0;
    int ncHadrons=0;
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
             //std::cout << "pdg ids = " << jet.constituents().at(j).user_index() << std::endl;
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
    for (unsigned int j = 0; j < HFparticles.size(); j++){
	if(deltaR(HFparticles[j].eta(),HFparticles[j].phi(),curjet.eta(),curjet.phi())<(RPARAM-0.1*RPARAM)){
		//std::cout<<nbHadrons<< "   "<<ncHadrons<<"  "<<deltaR(HFparticles[j].eta(),HFparticles[j].phi(),curjet.eta(),curjet.phi())<<"  "<<HFparticles[j].user_index()<<std::endl;
		if(abs(HFparticles[j].user_index())/100==4 || (abs(HFparticles[j].user_index())-10000)/100==4 || (abs(HFparticles[j].user_index())-20000)/100==4) ncHadrons++;
		else if (abs(HFparticles[j].user_index())/100==5 || (abs(HFparticles[j].user_index())-10000)/100==5 || (abs(HFparticles[j].user_index())-20000)/100==5 || abs(HFparticles[j].user_index())/1000==5) nbHadrons++;
//		else std::cout<< HFparticles[j].user_index()<<std::endl;
	    }	
	
	}
    //std::cout<<nbHadrons<< "   "<<ncHadrons<<std::endl;
    j_nbHadrons.push_back(nbHadrons);
    j_ncHadrons.push_back(ncHadrons);
    // groomers/taggers
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    
    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0;
    fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdb2(2.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdm1(-1.0, zcut_sd, mu_sd);
    fastjet::RecursiveSoftDrop soft_drop_rsd(0.0, zcut_sd, RPARAM);
    
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
    /*
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b1 (1.0,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b2 (1.0,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    */

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
	double beta=1;
	std::vector<fastjet::PseudoJet> lConst = curjet.constituents();
	fECF->calcECFN(beta,lConst,true);
        // energy correlator     
        j_zlogz.push_back( zlogz );   
        j_c1_b0.push_back( ECF_C1_b0(curjet) );
        j_c1_b1.push_back( ECF_C1_b1(curjet) );
        j_c1_b2.push_back( ECF_C1_b2(curjet) );
        j_c2_b1.push_back( ECF_C2_b1(curjet) );
        j_c2_b2.push_back( ECF_C2_b2(curjet) );
	j_d2_b1.push_back( ECF_D2_b1(curjet) );
        j_d2_b2.push_back( ECF_D2_b2(curjet) );

	j_d2_a1_b1.push_back( fECF->D2()        );
	j_m2_b1   .push_back( fECF->M2()        );
        j_n2_b1   .push_back( fECF->N2()        );
	j_m3_b1   .push_back( fECF->M3()        );
        j_n3_b1   .push_back( fECF->N3()        );

	beta=2;
	lConst = curjet.constituents();
	fECF->calcECFN(beta,lConst);
	j_d2_a1_b2.push_back( fECF->D2()        );
        j_m2_b2.push_back   ( fECF->M2()        );
        j_n2_b2.push_back   ( fECF->M2()        );
	j_m3_b2.push_back   ( fECF->M3()        );
        j_n3_b2.push_back   ( fECF->N3()        );
	/*
	pAddJet->e2_b1      = float(fECF->manager->ecfns["2_2"]);
	pAddJet->e3_b1      = float(fECF->manager->ecfns["3_3"]);
	pAddJet->e3_v1_b1   = float(fECF->manager->ecfns["3_1"]);
	pAddJet->e3_v2_b1   = float(fECF->manager->ecfns["3_2"]);
	pAddJet->e4_v1_b1   = float(fECF->manager->ecfns["4_1"]);
	pAddJet->e4_v2_b1   = float(fECF->manager->ecfns["4_2"]);
	*/
	
        // Groomed variables
        fastjet::PseudoJet mmdtjet=soft_drop_mmdt(curjet);
        fastjet::PseudoJet rsdjet=soft_drop_rsd(curjet);

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

	beta=1;
	lConst = mmdtjet.constituents();
	fECF->calcECFN(beta,lConst,true);
	j_d2_a1_b1_mmdt.push_back( fECF->D2()        );
	j_m2_b1_mmdt   .push_back( fECF->M2()        );
        j_n2_b1_mmdt   .push_back( fECF->N2()        );
	j_m3_b1_mmdt   .push_back( fECF->M3()        );
        j_n3_b1_mmdt   .push_back( fECF->N3()        );

	beta=2;
	lConst = mmdtjet.constituents();
	fECF->calcECFN(beta,lConst);
	j_d2_a1_b2_mmdt.push_back( fECF->D2()        );
        j_m2_b2_mmdt.push_back   ( fECF->M2()        );
        j_n2_b2_mmdt.push_back   ( fECF->M2()        );
	j_m3_b2_mmdt.push_back   ( fECF->M3()        );
        j_n3_b2_mmdt.push_back   ( fECF->N3()        );

        j_mass_trim.push_back( trimmer1( curjet ).m() );
        j_mass_prun.push_back( pruner1( curjet ).m() );    
        j_mass_mmdt.push_back( mmdtjet.m() );
        j_mass_rsd.push_back( rsdjet.m() );
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
        j_m3_b1.push_back( -99 );
        j_m3_b2.push_back( -99 );
        j_n3_b1.push_back( -99 );
        j_n3_b2.push_back( -99 );

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
        j_m3_b1_mmdt.push_back( -99 );
        j_m3_b2_mmdt.push_back( -99 );
        j_n3_b1_mmdt.push_back( -99 );
        j_n3_b2_mmdt.push_back( -99 );
        
        j_mass_trim.push_back( -99 );
        j_mass_prun.push_back( -99 );
        j_mass_mmdt.push_back( -99 );
        j_mass_rsd.push_back( -99 );
        j_mass_sdb2.push_back( -99 );
        j_mass_sdm1.push_back( -99 );
        
        j_multiplicity.push_back( -99 );

    }
}

// ----------------------------------------------------------------------------------
void smearJetPt(fastjet::PseudoJet &jet, Double_t shiftChargedScale, Double_t shiftPhotonScale, Double_t shiftHadronScale) {

    Double_t energyResolution=0;
    Double_t energyScale=1.;
    for(auto& it: ids) { 
        std::vector<int> tidSet = it.second;

        if (std::find(tidSet.begin(), tidSet.end(), jet.user_index()) != tidSet.end()) {
            if(it.first=="c") {                                    // -------------- charged hadron resolution
                //energyResolution=sqrt(pow(0.0001*jet.pt(),2)+pow(0.005,2)); // CMS TDR
                energyResolution=sqrt(pow(0.00025*jet.pt(),2)+pow(0.015,2)); // Delphes CMS tuning https://github.com/sethzenz/Delphes/blob/master/Cards/CMS_Phase_I_NoPileUp.tcl#L159
                if (shiftChargedScale>0) energyScale=shiftChargedScale;
            } else if (it.first=="p") {                              // ---------------------- photon resolution
                //energyResolution=sqrt(pow(0.027/sqrt(jet.e()),2)+pow(0.15/jet.e(),2)+pow(0.005,2)); // CMS TDR
                //energyResolution=sqrt(pow(0.027/sqrt(1.6*jet.e()),2)+pow(0.15/1.6/jet.e(),2)+pow(1.6*0.005,2)); // CMS TDR + factor 1.6=1+0.6 for the fact that 60% charged hadrons (leaving 30% in ECAL) are substracted before 30% photons are reconstructed
                energyResolution=sqrt(pow(0.027/sqrt(jet.e()/1.6),2)+pow(0.15/jet.e()*1.6,2)+pow(0.005,2)); // CMS TDR + factor 1.6=1+0.6 for the fact that 60% charged hadrons (leaving 30% in ECAL) are substracted before 30% photons are reconstructed
                //energyResolution=sqrt(pow(0.042/jet.e(),2)+pow(0.1/sqrt(jet.e()),2)+pow(0.005,2)); // CMS tuning https://github.com/cms-met/cmssw/blob/b742cc16aff1915b275cf0847dcff93aa6deab14/RecoMET/METProducers/python/METSigParams_cfi.py#L37
                if (shiftPhotonScale>0) energyScale=shiftPhotonScale;
            } else if (it.first=="n") {                              // -------------- neutral hadron resolution
                //energyResolution=sqrt(pow(1.20/sqrt(jet.e()),2)+pow(0.05,2)); // CMS JET JINST
		//energyResolution=sqrt(pow(1.20/sqrt(7.*jet.e()),2)+pow(7.*0.05,2)); // CMS JET JINST + factor 7=1+6 for the fact that 60% charged hadrons are substracted before 10% neutral hadrons are reconstructed
                energyResolution=sqrt(pow(1.20/sqrt(jet.e()/7.),2)+pow(0.05,2)); // CMS JET JINST + factor 7=1+6 for the fact that 60% charged hadrons are substracted before 10% neutral hadrons are reconstructed
                //energyResolution=sqrt(pow(0.41/jet.e(),2)+pow(0.52/sqrt(jet.e()),2)+pow(0.25,2)); // CMS tuning https://github.com/cms-met/cmssw/blob/b7T42cc16aff1915b275cf0847dcff93aa6deab14/RecoMET/METProducers/python/METSigParams_cfi.py#L37
                if (shiftHadronScale>0) energyScale=shiftHadronScale;
            } else if (it.first=="l") {                              // ---------------------- lepton resolution
                //energyResolution=sqrt(pow(0.0001*jet.pt(),2)+pow(0.005,2)); // CMS TDR
                energyResolution=sqrt(pow(0.00025*jet.pt(),2)+pow(0.015,2)); // Delphes CMS tuning https://github.com/sethzenz/Delphes/blob/master/Cards/CMS_Phase_I_NoPileUp.tcl#L159
                if (shiftChargedScale>0) energyScale=shiftChargedScale;
            } else {                                                // ------------------------ everything else
                std::cout << "WARNING: Input jet has unlisted PDG ID!" << std:: endl;
                energyResolution=0; // assume perfect resolution for everything else (?)
            }
            // std::cout << it.first << " " << energyResolution << std::endl;
        }
    }

    float resFudgeFactor = 1.;
    bool nosmear = 0.;
    Double_t smearedPt = std::max(1e-10,smearDist->Gaus(1,energyResolution*resFudgeFactor));

    if (nosmear) smearedPt = 1.;
    smearedPt*=energyScale;

    jet.reset_momentum(jet.px() * smearedPt,
                       jet.py() * smearedPt,
                       jet.pz() * smearedPt,
                       jet.e()  * smearedPt);

}

// ----------------------------------------------------------------------------------
std::vector<fastjet::PseudoJet> discretizeEvent(std::vector<fastjet::PseudoJet> &particles, bool discretizeEcal, Double_t numberOfPileup, Double_t maxChargedPt, Double_t maxChargedDr, Double_t trackingEfficiency)
{

    static Double_t etaMin = -3,
                    etaMax = 3,
                    etaRes = 0.087/sqrt(12), // CMS TDR, factor 3 since particle flow cluster algorithm reaches this precision
                    phiMin = 0, 
                    phiMax = 2*TMath::Pi(),
                    phiRes = 0.087/sqrt(12); // CMS TDR, factor 3 since particle flow cluster algorithm reaches this precision
    static Int_t etaNBins = TMath::Floor((etaMax - etaMin)/etaRes), 
                 phiNBins = TMath::Floor((phiMax - phiMin)/phiRes);

    TH2D* hcalGrid = new TH2D("hcalGrid","hcalGrid",etaNBins,etaMin,etaMax,phiNBins,phiMin,phiMax);                   

    static Double_t etaResEcal = 0.017, // CMS ECAL JINST
                    phiResEcal = 0.017; // CMS ECAL JINST
    static Int_t etaNBinsEcal = TMath::Floor((etaMax - etaMin)/etaResEcal), // CMS ECAL JINST
                 phiNBinsEcal = TMath::Floor((phiMax - phiMin)/phiResEcal); // CMS ECAL JINST

    TH2D* ecalGrid = new TH2D("ecalGrid","ecalGrid",etaNBinsEcal,etaMin,etaMax,phiNBinsEcal,phiMin,phiMax);  
                 
    static Double_t pileupEnergyInCell = numberOfPileup * 0.3 *etaRes*phiRes; // neutral pileup density rho = 0.3 GeV / unit area / pileup interaction
    for (int i = 0; i < etaNBins; ++i)
        for (int j = 0; j < phiNBins; ++j)
            hcalGrid->SetBinContent(i+1,j+1,pileupEnergyInCell);

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
                if(track&&(maxChargedPt<10000)&&(particles[i].pt()>maxChargedPt))
                   track=false;
                if((track)||(it.first=="l")||((!discretizeEcal) && (it.first=="p")))
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
                    // std::cout << "n par id = " << particles[i].user_index() << "," << particles[i].e() << "," << ieta << "," << iphi << "," << cureta << "," << curphi << std::endl;
                }
                else if ((it.first=="n")||(!track)){
                    double curphi = particles[i].phi();
                    double cureta = particles[i].eta();
                    int ieta = hcalGrid->GetXaxis()->FindBin( cureta );
                    int iphi = hcalGrid->GetYaxis()->FindBin( curphi );
                    hcalGrid->SetBinContent( ieta,iphi, hcalGrid->GetBinContent(ieta,iphi)+particles[i].e() );
                    // std::cout << "n par id = " << particles[i].user_index() << "," << particles[i].e() << "," << ieta << "," << iphi << "," << cureta << "," << curphi << std::endl;
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
                // std::cout << "celle = " << celle << ", cellpt = " << cellpt << ", orig e = " << curbincontent << ", cell eta = " << celleta << std::endl;
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
                // std::cout << "celle = " << celle << ", cellpt = " << cellpt << ", orig e = " << curbincontent << ", cell eta = " << celleta << std::endl;
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

std::vector<fastjet::PseudoJet> discretizeJet(fastjet::PseudoJet jet, 
                                              bool clusterJet) 
{
    static fastjet::JetDefinition   jetDef(fastjet::antikt_algorithm, RPARAM);    
    static fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    static fastjet::AreaDefinition  fjAreaDefinition( fastjet::active_area, fjActiveArea );
    static Double_t etaMin = -10,
                    etaMax = 10,
                    etaRes = 0.005,
                    phiMin = 0, 
                    phiMax = 2*TMath::Pi(),
                    phiRes = 0.005; 
    static Int_t etaNBins = TMath::Floor((etaMax - etaMin)/etaRes), 
                 phiNBins = TMath::Floor((phiMax - phiMin)/phiRes);
    
    std::map<TString,std::pair<TH2D*,TH2D*>> resolutionHistos =
    { { "n" , std::pair<TH2D*,TH2D*>(new TH2D("","",etaNBins,etaMin,etaMax,
                                                    phiNBins,phiMin,phiMax),
                                     new TH2D("","",etaNBins,etaMin,etaMax,
                                                    phiNBins,phiMin,phiMax)) } };


    std::vector< fastjet::PseudoJet > newPar;
    
    // collect constituents
    for(int iConst=0; iConst < jet.constituents().size(); iConst++) {
        fastjet::PseudoJet constituent = jet.constituents().at(iConst);

        bool isDiscretizedType = false;
        for(auto const& it: ids) {
            if (resolutionHistos[it.first].first == 0) continue;
            std::vector<int> tidSet = it.second;

            if (std::find(tidSet.begin(), tidSet.end(), constituent.user_index()) != tidSet.end()) {
               isDiscretizedType=true;
               resolutionHistos[it.first].first->Fill(constituent.eta(),
                                                      constituent.phi(),
                                                      constituent.e());
               resolutionHistos[it.first].second->Fill(constituent.eta(),
                                                       constituent.phi(),
                                                        pow(constituent.px(),2)
                                                       +pow(constituent.py(),2)
                                                       +pow(constituent.pz(),2));
            }
        }

        if(!isDiscretizedType) {
            newPar.push_back(constituent);
        }
    }

    // build new constituents
    for(auto& it: resolutionHistos) {
        if (it.second.first == 0) continue;
        TH2D* EHisto = it.second.first;
        TH2D* pTHisto= it.second.second;

        // loop over bins
        for(int iEta=1; iEta <= etaNBins; iEta++) {
        for(int iPhi=1; iPhi <= phiNBins; iPhi++) {
            Int_t binNum    = EHisto->GetBin(iEta,iPhi);
            Double_t binEta = EHisto->GetXaxis()->GetBinCenter(iEta);
            Double_t binPhi = EHisto->GetYaxis()->GetBinCenter(iPhi);
            Double_t binE   = EHisto->GetBinContent(binNum);
            Double_t binp   = sqrt(pTHisto->GetBinContent(binNum));
            if(binE == 0 && binp == 0) continue;

            Double_t px = binp*TMath::Cos(binPhi),
                     py = binp*TMath::Sin(binPhi),
                     pz = binp*sinh(binEta);

            // reconstruct PseudoJet using pE 4-vector
            // (this is the only way in fastjet unless we assume m=0)
            fastjet::PseudoJet newConst(px,py,pz,binE);
            newConst.set_user_index(ids[it.first][0]); // just steal some id for now
            newPar.push_back(newConst);
        }}

        delete EHisto;
        delete pTHisto;
    }
    
    for(auto& it: resolutionHistos) {
        delete it.second.first;
        delete it.second.second;
    }

    // build jet or just return constituents
    if(clusterJet) {
        if(newPar.size() > 0) {
            fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newPar, jetDef, fjAreaDefinition);
            std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));
            if(out_jets.size()>0) {
              thisClustering->delete_self_when_unused();
            } else delete thisClustering;
            return out_jets;
        } else {
            fastjet::PseudoJet nulljet = fastjet::PseudoJet(0,0,0,0);
            newPar.push_back(nulljet);
            return newPar;
        }
    } else {
        return newPar;
    }
}
