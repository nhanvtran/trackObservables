# ----------------------------------------
# inputs 

FASTJETLOCATION=/uscms_data/d2/ntran/physics/ttbarMela/fastjet/fastjet-install
#/afs/cern.ch/user/h/hinzmann/workspace/Wmass/trackObservables/processing/fastjet/fastjet-install
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

HEPMCLOCATION=/uscms_data/d2/ntran/physics/ttbarMela/HepMC/install
#/afs/cern.ch/user/h/hinzmann/workspace/Wmass/trackObservables/processing/HepMC/HepMC-install
HEPMCINC :=-I$(HEPMCLOCATION)/include
HEPMCLIB :=-L$(HEPMCLOCATION)/lib  -lHepMC

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

all: anaSubstructure

anaSubstructure: anaSubstructure.cpp
#	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lRecursiveTools -lEnergyCorrelator anaSubstructure.cpp 
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) $(HEPMCINC) anaSubstructure.cpp $(HEPMCLIB) -lEnergyCorrelator -lRecursiveTools -lNsubjettiness -lz puppiContainer.cc PU14.cc EventMixer.cc CmdLine.cc ParticleIDMethods.cc ParticleName.cc EventSource.cc zfstream.cc Version.cc EnergyCorrelations.cc RecursiveSoftDrop.cc

clean:
	rm anaSubstructure
