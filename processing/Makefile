# ----------------------------------------
# inputs 

FASTJETLOCATION=./fastjet-install
FASTJETINC :=`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`

ROOTINC := `root-config --glibs --cflags`

# ----------------------------------------

all: anaSubstructure

anaSubstructure: anaSubstructure.cpp
#	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) -lRecursiveTools -lEnergyCorrelator anaSubstructure.cpp 
	c++ -o anaSubstructure $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) anaSubstructure.cpp -lEnergyCorrelator -lRecursiveTools -lNsubjettiness

clean:
	rm anaSubstructure