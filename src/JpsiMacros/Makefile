ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  
LD           = c++
LDFLAGS       = -O2 -g
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
ifeq ($(MACOSX_MINOR),4)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.4 c++
else
ifeq ($(MACOSX_MINOR),3)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.3 c++
else
UNDEFOPT      = suppress
LD            = c++
endif
endif

endif

CXXFLAGS      += $(ROOTCFLAGS)

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)
OUTLIB	         = ../lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../lib/

$(OUTLIB)JPsiTreeBase.o: $(INCLUDEDIR)/JPsiTreeBase.C $(INCLUDEDIR)/finalJPsiAnalysis.cc 
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)JPsiTreeBase.o $<

$(OUTLIB)finalJPsiAnalysis.o: $(INCLUDEDIR)/finalJPsiAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)finalJPsiAnalysis.o $<

$(OUTLIB)MakeDataSet.o: $(INCLUDEDIR)/MakeDataSet.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)MakeDataSet.o $<

$(OUTLIB)RooHistPdfConv.o: $(INCLUDEDIR)/RooHistPdfConv.cxx
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/libRooHistPdfConv.o $<

JPsiApp:  $(INCLUDEDIR)JPsiApp.C $(OUTLIB)JPsiTreeBase.o  
	$(CXX) $(CXXFLAGS) -o JPsiApp $(OUTLIB)/*.o $(GLIBS) $ $<

JPsiFitApp:  $(INCLUDEDIR)JPsiFitApp.C $(OUTLIB)JPsiTreeBase.o  
	$(CXX) $(CXXFLAGS) -o JPsiFitApp $(OUTLIB)/*.o $(GLIBS) $ $<

FitMass:  $(INCLUDEDIR)fitMass.cxx
	$(CXX) $(CXXFLAGS) -o FitMass $(OUTLIB)/*.o $(GLIBS) $ $<

FitMassSim:  $(INCLUDEDIR)fitMass_simul.cxx
	$(CXX) $(CXXFLAGS) -o FitMassSim $(OUTLIB)/*.o $(GLIBS) $ $<

ShowTruth:  $(INCLUDEDIR)showTruth.cxx
	$(CXX) $(CXXFLAGS) -o ShowTruth $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2D:    $(INCLUDEDIR)fit2D.cxx $(OUTLIB)RooHistPdfConv.o  
	$(CXX) $(CXXFLAGS) -o Fit2D $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DGT:    $(INCLUDEDIR)fit2DGT.cxx
	$(CXX) $(CXXFLAGS) -o Fit2DGT $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DSim: $(INCLUDEDIR)fit2D_sim.cxx
	$(CXX) $(CXXFLAGS) -o Fit2DSim $(OUTLIB)/*.o $(GLIBS) $ $<

FitSmearFactor: $(INCLUDEDIR)fitSmearFactor.cxx
	$(CXX) $(CXXFLAGS) -o FitSmearFactor $(OUTLIB)/*.o $(GLIBS) $ $<

FitSidebands:  $(INCLUDEDIR)fitSidebands.cxx
	$(CXX) $(CXXFLAGS) -o FitSidebands $(OUTLIB)/*.o $(GLIBS) $ $<

FitLifetime:  $(INCLUDEDIR)fitLifetime.cxx
	$(CXX) $(CXXFLAGS) -o FitLifetime $(OUTLIB)/*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f $(OUTLIB)*.so
	rm -f JPsiApp JPsiFitApp

all:  JPsiApp JPsiFitApp Fit2D FitMass ShowTruth FitMassSim
