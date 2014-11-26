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
     //isVBF( pset.getUntrackedParameter<bool>("isVBF",false) ),
     //lookFor( pset.getUntrackedParameter<int>("lookForGenProcID",0) ),
     fOutputFile(0)
{
}

void LHEPythiaEventAnalyzer::beginJob()
{
 
   nevent = 0;
   neventpt0 = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
 
   pTZ = new TH1D("pTZ", "pT of Z boson", 21, 0., 140.) ;
   nJets = new TH1D("nJets", "# of GenJets with pT > 20", 5, -0.5, 4.5) ; 
   pTjet1 = new TH1D("pTjet1", "pT of leading GenJet", 18, 10., 100.) ;
   pTjet2 = new TH1D("pTjet2", "pT of second leading GenJet", 18, 10., 100.) ;
   pTjet3 = new TH1D("pTjet3", "pT of third leading GenJet", 18, 0., 100.) ;
 
   return ;
}
 
void LHEPythiaEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   e.getByLabel(std::string("genParticles"), genParticles);
   TLorentzVector mu[2];
   int iMu = 0;

   for(std::vector<reco::GenParticle>::const_iterator genParticle=genParticles->begin(); genParticle!=genParticles->end(); ++genParticle){
     if( abs(genParticle->pdgId()==13) && genParticle->status()==1 ) {
       mu[iMu].SetPxPyPzE(genParticle->px(),genParticle->py(),
                          genParticle->pz(),genParticle->energy());
       iMu++;
       if (iMu == 2) break;
     }
   }
   TLorentzVector mumu = mu[0] + mu[1];
   pTZ->Fill(mumu.Perp());

   edm::Handle<reco::GenJetCollection> jetsgen;
   e.getByLabel("ak5GenJets", jetsgen);
   // const reco::GenJetCollection genJets = jetsgen.product();
   reco::GenJetCollection* genJetsSorted = new reco::GenJetCollection();

   for(reco::GenJetCollection::const_iterator genJet=jetsgen->begin(); genJet!=jetsgen->end(); ++genJet){
     if( genJet->pt() > 20. ) {     
       TLorentzVector jet4;
       jet4.SetPxPyPzE(genJet->px(),genJet->py(),
		       genJet->pz(),genJet->energy());
       if (jet4.DeltaR(mu[0]) > 0.5 && jet4.DeltaR(mu[1]) > 0.5)  
	 genJetsSorted->push_back(*genJet); 
     }
   }
    
   std::sort(genJetsSorted->begin(), genJetsSorted->end(), pTComparator_);
  
   nJets->Fill(genJetsSorted->size());
   if (genJetsSorted->size() > 0) pTjet1->Fill(genJetsSorted->at(0).pt());
   if (genJetsSorted->size() > 1) pTjet2->Fill(genJetsSorted->at(1).pt());
   if (genJetsSorted->size() > 2) pTjet3->Fill(genJetsSorted->at(2).pt()); 

   return ;   
}

void LHEPythiaEventAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(pTZ);	   
  Hlist.Add(nJets) ;
  Hlist.Add(pTjet1);
  Hlist.Add(pTjet2); 
  Hlist.Add(pTjet3) ;
  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n"; 
  cout << "N_events(pT = 0) = " << neventpt0 << "\n"; 
  return ;
}
 
DEFINE_FWK_MODULE(LHEPythiaEventAnalyzer);
