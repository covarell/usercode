//#include "TROOT.h"
//#include "TSystem.h"
//#include "MELA.C"

void Run_MELA(char* inputFile, float minMzz = 100., float maxMzz = 1000., bool withPT = false, bool withY = false)
{
  //gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
  gROOT->ProcessLine(".L ../PDFs/RooXZsZs_5D.cxx+");                
  gROOT->ProcessLine(".L ../src/AngularPdfFactory.cc+");                 
  gROOT->ProcessLine(".L ../PDFs/RooqqZZ_JHU.cxx+");                         
  gROOT->ProcessLine(".L ../PDFs/RooTsallis.cc+");
  gROOT->ProcessLine(".L ../PDFs/RooRapidityBkg.cxx+");
  gROOT->ProcessLine(".L ../PDFs/RooRapiditySig.cxx+");
  gROOT->ProcessLine(".L MELA.C+");
  addDtoTree(inputFile, minMzz, maxMzz, withPT, withY);

}
