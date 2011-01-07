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
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/cms/cmssw/CMSSW_3_6_3/external/slc5_ia32_gcc434/lib -lMinuit -lRooFit
# NGLIBS        += -lMinuit -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)
# Hack here, because RooFit is removed from ROOT:
CXX	         += -I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms6/include
OUTLIB	         = ../../../../../lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../../../../../lib/

## $(OUTLIB)JPsiTreeBase.o: $(INCLUDEDIR)/JPsiTreeBase.C $(INCLUDEDIR)/finalJPsiAnalysis.cc 
##	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)JPsiTreeBase.o $<

## $(OUTLIB)finalJPsiAnalysis.o: $(INCLUDEDIR)/finalJPsiAnalysis.cc
##	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)finalJPsiAnalysis.o $<

## $(OUTLIB)MakeDataSet.o: $(INCLUDEDIR)/MakeDataSet.cc
##	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)MakeDataSet.o $<

RooHistPdfConv.o: $(INCLUDEDIR)/RooHistPdfConv.cxx
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/libRooHistPdfConv.o $<

## JPsiApp:  $(INCLUDEDIR)JPsiApp.C $(OUTLIB)JPsiTreeBase.o  
##	$(CXX) $(CXXFLAGS) -o JPsiApp $(OUTLIB)/*.o $(GLIBS) $ $<

TreeToDat:  $(INCLUDEDIR)tree_to_datasets.cxx  
	$(CXX) $(CXXFLAGS) -o TreeToDat $(GLIBS) $ $<

FitMass:  $(INCLUDEDIR)fitMass.cxx
	$(CXX) $(CXXFLAGS) -o FitMass $(GLIBS) $ $<

FitMassData:  $(INCLUDEDIR)fitMassData.cxx
	$(CXX) $(CXXFLAGS) -o FitMassData $(GLIBS) $ $<

FitMassDataRange:  $(INCLUDEDIR)fitMassData_ranges.cxx
	$(CXX) $(CXXFLAGS) -o FitMassDataRange $(GLIBS) $ $<

FitMassData2S:  $(INCLUDEDIR)fitMassDataPsi2S.cxx
	$(CXX) $(CXXFLAGS) -o FitMassData2S $(GLIBS) $ $<

FitMassSmall:  $(INCLUDEDIR)fitMassSmall.cxx
	$(CXX) $(CXXFLAGS) -o FitMassSmall $(GLIBS) $ $< 

FitMassRange:  $(INCLUDEDIR)fitMass_ranges.cxx
	$(CXX) $(CXXFLAGS) -o FitMassRange $(GLIBS) $ $<

FitMassSim:  $(INCLUDEDIR)fitMass_simul.cxx
	$(CXX) $(CXXFLAGS) -o FitMassSim $(GLIBS) $ $<

Fit2D:    $(INCLUDEDIR)fit2D.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2D $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DData:    $(INCLUDEDIR)fit2DData.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2DData $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DDataRange:    $(INCLUDEDIR)fit2DData_ranges.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2DDataRange $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DDataPsiP:    $(INCLUDEDIR)fit2DData_psip.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2DDataPsiP $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DDataSyst:    $(INCLUDEDIR)fit2DData_syst.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2DDataSyst $(OUTLIB)/*.o $(GLIBS) $ $<

Toy2DDataRange:    $(INCLUDEDIR)toy2DData_ranges.cxx   
	$(CXX) $(CXXFLAGS) -o Toy2DDataRange $(OUTLIB)/*.o $(GLIBS) $ $<

Fit2DRange:    $(INCLUDEDIR)fit2D_ranges.cxx   
	$(CXX) $(CXXFLAGS) -o Fit2DRange $(OUTLIB)/*.o $(GLIBS) $ $<

ShowTruth:  $(INCLUDEDIR)showTruth.cxx
	$(CXX) $(CXXFLAGS) -o ShowTruth $(OUTLIB)/*.o $(GLIBS) $ $<

## Fit2DGT:    $(INCLUDEDIR)fit2DGT.cxx
##	$(CXX) $(CXXFLAGS) -o Fit2DGT $(OUTLIB)/*.o $(GLIBS) $ $<

## Fit2DSim: $(INCLUDEDIR)fit2D_sim.cxx
##	$(CXX) $(CXXFLAGS) -o Fit2DSim $(OUTLIB)/*.o $(GLIBS) $ $<

## FitSmearFactor: $(INCLUDEDIR)fitSmearFactor.cxx
##	$(CXX) $(CXXFLAGS) -o FitSmearFactor $(OUTLIB)/*.o $(GLIBS) $ $<

FindMeans:  $(INCLUDEDIR)findMeans.cxx
	$(CXX) $(CXXFLAGS) -o FindMeans $(GLIBS) $ $<

Compare:  $(INCLUDEDIR)compareWithMC.cxx
	$(CXX) $(CXXFLAGS) -o Compare $(GLIBS) $ $<

FitLifetimePrompt:  $(INCLUDEDIR)fitLifetimePrompt.cxx
	$(CXX) $(CXXFLAGS) -o FitLifetimePrompt $(GLIBS) $ $<

FitLifetimeNonPrompt:  $(INCLUDEDIR)fitLifetimeNonPrompt.cxx
	$(CXX) $(CXXFLAGS) -o FitLifetimeNonPrompt $(OUTLIB)/*.o $(GLIBS) $ $<

FitLifetimeBkg:  $(INCLUDEDIR)fitLifetimeBkg.cxx
	$(CXX) $(CXXFLAGS) -o FitLifetimeBkg $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f $(OUTLIB)*.so
##	rm -f JPsiApp JPsiFitApp

all:  Fit2D FitMass FitMassSim
