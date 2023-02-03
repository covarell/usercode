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

#include "GenMacros/GenMacros/test/LHEEventAnalyzer.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
LHEEventAnalyzer::LHEEventAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     whichWeight( pset.getUntrackedParameter<int>("whichWeight",-1)),
     xsection( pset.getUntrackedParameter<double>("xSection",1.)),
     fOutputFile(0)
{

   nevent = 0;
   lhep_token = consumes< LHEEventProduct >(InputTag(theSrc));
   
   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   fourlMass  = new TH1D(  "fourlMass", "Four-lepton mass", 100,  70., 970. ) ;  

   //mult = new TH1D( "mult", "Parton multiplicity" , 5, -0.5 , 4.5 ) ; 
   fourlMass->Sumw2();
}

LHEEventAnalyzer::~LHEEventAnalyzer()
{
  TObjArray Hlist(0);
  fourlMass->Scale(1./float(nevent));
  Hlist.Add(fourlMass);	   
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
   float mass = 0.;
   TLorentzVector fourl(0.,0.,0.,0.);

   Handle< LHEEventProduct > EvtHandle ;
   e.getByToken( lhep_token , EvtHandle ) ;

   float weight = xsection;
   if (whichWeight >= 0) weight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 
   if (nevent < 20) cout << "weight = " << weight << endl;

   for (int i = 0; i < EvtHandle->hepeup().NUP; ++i) {
     if (EvtHandle->hepeup().ISTUP[i] != 1) {
           continue;
     }
     if (abs(EvtHandle->hepeup().IDUP[i]) == 11 || abs(EvtHandle->hepeup().IDUP[i]) == 13) {
       TLorentzVector temp; 
       temp.SetPxPyPzE(EvtHandle->hepeup().PUP[i][0],EvtHandle->hepeup().PUP[i][1],
		       EvtHandle->hepeup().PUP[i][2],EvtHandle->hepeup().PUP[i][3]);
       fourl += temp;
     }

   }
   fourlMass->Fill(fourl.M(),weight);

   return ;   
}

void LHEEventAnalyzer::endJob()
{
}
 
DEFINE_FWK_MODULE(LHEEventAnalyzer);
