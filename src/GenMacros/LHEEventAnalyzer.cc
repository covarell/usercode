// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      LHEEventAnalyzer
// 
//
// Description: Module to analyze Pythia-EvtGen HepMCproducts
//
//
// Original Author:  Roberto Covarelli
//         Created:  April 26, 2007
//

#include <iostream>
#include <fstream>
 
// #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObjArray.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Configuration/GenProduction/test/LHEEventAnalyzer.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
LHEEventAnalyzer::LHEEventAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     whichWeight( pset.getUntrackedParameter<int>("whichWeight",-1)),
     fOutputFile(0)
{

   nevent = 0;
 
   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  

   mult = new TH1D( "mult", "Parton multiplicity" , 5, -0.5 , 4.5 ) ; 
   mult->Sumw2();
}

LHEEventAnalyzer::~LHEEventAnalyzer()
{
  TObjArray Hlist(0);
  Hlist.Add(mult);	   
  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n";  
}

void LHEEventAnalyzer::beginJob()
{
}
 
void LHEEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;
   int npart = 0;

   Handle< LHEEventProduct > EvtHandle ;
   e.getByLabel( theSrc , EvtHandle ) ;

   float weight = EvtHandle->hepeup().XWGTUP;
   if (whichWeight >= 0) weight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 

   for (int i = 0; i < EvtHandle->hepeup().NUP; ++i) {
     if (EvtHandle->hepeup().ISTUP[i] != 1) {
       //cout << reader->hepeup.ISTUP[i] << ", " << reader->hepeup.IDUP[i] << endl;
       continue;
     }
     // std::cout << "PDG id LHE = " << EvtHandle->hepeup().IDUP[i] << std::endl;
     if (EvtHandle->hepeup().IDUP[i] == 21 || abs(EvtHandle->hepeup().IDUP[i]) < 6) npart++;

     mult->Fill(npart,weight);
   }
   
   return ;   
}

void LHEEventAnalyzer::endJob()
{
}
 
DEFINE_FWK_MODULE(LHEEventAnalyzer);
