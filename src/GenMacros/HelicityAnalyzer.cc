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

#include "Configuration/GenProd/test/LHEPythiaEventAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


using namespace edm;
using namespace std;
// using namespace HepMC;
 
LHEPythiaEventAnalyzer::LHEPythiaEventAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     whichWeight( pset.getUntrackedParameter<int>("whichWeight",-1)),
     //isVBF( pset.getUntrackedParameter<bool>("isVBF",false) ),
     //lookFor( pset.getUntrackedParameter<int>("lookForGenProcID",0) ),
     fOutputFile(0)
{

  lhep_token = consumes< LHEEventProduct >(InputTag(theSrc));
  genp_token = consumes<std::vector<reco::GenParticle> >(InputTag("prunedGenParticles"));
}

void LHEPythiaEventAnalyzer::beginJob()
{
 
   nevent = 0;
   neventpt0 = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
   helAngleKPlus = new TH1D("helAngleKPlus", "cosine of K+ helicity angle", 80, -1., 1.) ;
   // helAngleKPlus->Sumw2(); 

   return ;
}
 
void LHEPythiaEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   e.getByToken(genp_token, genParticles);
   TLorentzVector mu[2];

   for(std::vector<reco::GenParticle>::const_iterator genParticle=genParticles->end(); genParticle!=genParticles->begin(); --genParticle){
     // std::cout << genParticle->pdgId() << std::endl;
     // last in the collection
     if( genParticle->pdgId()==321 && genParticle->mother() && genParticle->mother()->pdgId()==333) {
       const reco::GenParticle* phi = dynamic_cast<const reco::GenParticle*>(genParticle->mother());
       if( phi->mother() && phi->mother()->pdgId()==25) {
	  mu[1].SetPxPyPzE(genParticle->px(),genParticle->py(),
			   genParticle->pz(),genParticle->energy());
	  mu[0].SetPxPyPzE(phi->px(),phi->py(),
			   phi->pz(),phi->energy());
       }
     }
   }

   TVector3 phiBoost = mu[0].BoostVector();
   mu[1].Boost(-phiBoost);
   float theAngle = cos(mu[0].Vect().Angle(mu[1].Vect()));
   
   helAngleKPlus->Fill(theAngle);
   return ;   
}

void LHEPythiaEventAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(helAngleKPlus);
  Hlist.Write();
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n"; 
  cout << "N_events(pT = 0) = " << neventpt0 << "\n"; 
  return ;
}
 
DEFINE_FWK_MODULE(LHEPythiaEventAnalyzer);
