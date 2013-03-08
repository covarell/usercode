// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      LHEPythiaEventAnalyzer
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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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

#include "GeneratorInterface/ExternalDecays/test/LHEPythiaEventAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
LHEPythiaEventAnalyzer::LHEPythiaEventAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     isVBF( pset.getUntrackedParameter<bool>("isVBF",false) ),
     lookFor( pset.getUntrackedParameter<int>("lookForGenProcID",0) ),
     fOutputFile(0)
{
}

void LHEPythiaEventAnalyzer::beginJob()
{
 
   nevent = 0;
   neventpt0 = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
   string namept = "Pt_sig";
   string namey = "Y_sig";
   if (isVBF) {
     namept += "VBF";  
     namey += "VBF";
   }

   Pt_sigVBF = new TH2D( namept.c_str(), namept.c_str(), 420, 79.5, 499.5, 500, 0., 500.) ;
   Y_sigVBF = new TH2D( namey.c_str(), namey.c_str(), 420, 79.5, 499.5, 50, -5., 5.) ;
   Pt_bkg = new TH2D( "Pt_bkg", "Pt_bkg" , 420, 79.5, 499.5, 500, 0., 500.) ;
   PtOverM_bkg = new TH2D( "PtOverM_bkg", "Pt_bkg" , 420, 79.5, 499.5, 500, 0., 4.) ;
   Y_bkg = new TH2D( "Y_bkg" , "Y_bkg" , 420, 79.5, 499.5, 50, -5., 5.) ;
   return ;
}
 
void LHEPythiaEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   e.getByLabel(std::string("genParticles"), genParticles);
   // reco::GenParticle genHiggs;

   if (lookFor) {
     edm::Handle<GenEventInfoProduct> gen;
     e.getByLabel( "generator", gen );
     int processId = (int)(gen->signalProcessID());
     if ( abs(processId - lookFor) > 1 ) return;  // allow +/-1
   }

   for(std::vector<reco::GenParticle>::const_iterator genParticle=genParticles->begin(); genParticle!=genParticles->end(); ++genParticle){
     if(genParticle->pdgId()==25) {
       Pt_sigVBF->Fill(genParticle->mass(),genParticle->pt());
       Y_sigVBF->Fill(genParticle->mass(),genParticle->rapidity());
       break;
     }
   }
   return ;   
}

void LHEPythiaEventAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(Pt_sigVBF);	   
  Hlist.Add(Y_sigVBF) ;
  Hlist.Add(Pt_bkg);
  Hlist.Add(PtOverM_bkg); 
  Hlist.Add(Y_bkg) ;
  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n"; 
  cout << "N_events(pT = 0) = " << neventpt0 << "\n"; 
  return ;
}
 
DEFINE_FWK_MODULE(LHEPythiaEventAnalyzer);
