# ----------------------------------------
# inputs 

FASTJETLOCATION=/Users/ntran/Documents/Research/Ext/VHEPP/VHEPPStudies/fastjet/fastjet-install
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

all: anaSubstructure

anaSubstructure: anaSubstructure.cpp
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lRecursiveTools -lEnergyCorrelator anaSubstructure.cpp 

clean:
	rm anaSubstructure