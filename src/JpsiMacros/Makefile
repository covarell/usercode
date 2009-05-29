ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
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

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lRooFit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = ../lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ../lib/

$(OUTLIB)JPsiTreeBase.o: $(INCLUDEDIR)/JPsiTreeBase.C $(INCLUDEDIR)/finalJPsiAnalysis.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JPsiTreeBase.o $<

$(OUTLIB)finalJPsiAnalysis.o: $(INCLUDEDIR)/finalJPsiAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)finalJPsiAnalysis.o $<
$(OUTLIB)fitBfrac.o: $(INCLUDEDIR)/fitBfrac.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)fitBfrac.o $<

JPsiApp:  $(INCLUDEDIR)JPsiApp.C $(OUTLIB)JPsiTreeBase.o  
	$(CXX) $(CXXFLAGS) -o JPsiApp $(OUTLIB)/*.o $(GLIBS) $ $<

JPsiApp.clean:
	rm -f JPsiApp

JPsiFitApp:  $(INCLUDEDIR)JPsiFitApp.C $(OUTLIB)JPsiTreeBase.o  
	$(CXX) $(CXXFLAGS) -o JPsiFitApp $(OUTLIB)/*.o $(GLIBS) $ $<

JPsiFitApp.clean:
	rm -f JPsiFitApp

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f JPsiApp JPsiFitApp

all:  JPsiApp JPsiFitApp
