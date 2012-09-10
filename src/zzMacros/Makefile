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
#### CMSSW 4_2
# NGLIBS        += -L/afs/cern.ch/cms/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_4/external/slc5_amd64_gcc434/lib -lMinuit -lMinuit2 -lRooFit -lRooFitCore
#### CMSSW 5_2
NGLIBS        += -L/afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_2_5/external/slc5_amd64_gcc462/lib/ -lMinuit -lMinuit2 -lRooFit -lRooFitCore
# NGLIBS        += -lMinuit -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR)
# Hack here, because RooFit is removed from ROOT:
#### CMSSW 4_2
# CXX	         += -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include
#### CMSSW 5_2
CXX	         += -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms6/include
OUTLIB	         = lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

RooJetsPentaSpinTwo.o: $(INCLUDEDIR)RooJetsPentaSpinTwo.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/libRooJetsPentaSpinTwo.o $<

RooPentaSpinTwo.o: $(INCLUDEDIR)RooPentaSpinTwo.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/libRooPentaSpinTwo.o $<

FitAcc_ZZ:  $(INCLUDEDIR)fitAcc_ZZ.cxx
	$(CXX) $(CXXFLAGS) -o FitAcc_ZZ $(OUTLIB)/*.o $(GLIBS) $ $<

FitAccSimp_ZZ:  $(INCLUDEDIR)fitAcc_ZZ_simplif.cxx
	$(CXX) $(CXXFLAGS) -o FitAccSimp_ZZ $(OUTLIB)/*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f $(OUTLIB)*.so
##	rm -f JPsiApp JPsiFitApp

all: RooJetsPentaSpinTwo.o FitAcc_ZZ FitAccSimp_ZZ 
