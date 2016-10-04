# ----------------------------------------
# inputs 

FASTJETLOCATION=./fastjet/fastjet-install
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

all: anaSubstructure

anaSubstructure: anaSubstructure.cpp
#	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lRecursiveTools -lEnergyCorrelator anaSubstructure.cpp 
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) anaSubstructure.cpp -lEnergyCorrelator -lRecursiveTools -lNsubjettiness -lz puppiContainer.cc PU14.cc EventMixer.cc CmdLine.cc ParticleIDMethods.cc ParticleName.cc EventSource.cc zfstream.cc Version.cc

clean:
	rm anaSubstructure
