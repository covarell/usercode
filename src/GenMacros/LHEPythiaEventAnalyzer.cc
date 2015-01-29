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

#include "Configuration/GenProduction/test/LHEPythiaEventAnalyzer.h"
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
}

void LHEPythiaEventAnalyzer::beginJob()
{
 
   nevent = 0;
   neventpt0 = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
   mZ = new TH1D("mZ", "mass of Z boson", 21, 0., 140.) ;
   pTZ = new TH1D("pTZ", "pT of Z boson", 21, 0., 140.) ;
   nJets = new TH1D("nJets", "# of GenJets with pT > 20", 7, -0.5, 6.5) ; 
   pTjet1 = new TH1D("pTjet1", "pT of leading GenJet", 18, 10., 100.) ;
   pTjet2 = new TH1D("pTjet2", "pT of second leading GenJet", 18, 10., 100.) ;
   pTjet3 = new TH1D("pTjet3", "pT of third leading GenJet", 18, 10., 100.) ;
   		 
   pTZ_lowmass = new TH1D("pTZ_lowmass", "pT of Z boson", 21, 0., 140.) ;
   nJets_lowmass = new TH1D("nJets_lowmass", "# of GenJets with pT > 20", 5, -0.5, 4.5) ; 
   pTjet1_lowmass = new TH1D("pTjet1_lowmass", "pT of leading GenJet", 18, 10., 100.) ;
   pTjet2_lowmass = new TH1D("pTjet2_lowmass", "pT of second leading GenJet", 18, 10., 100.) ;
   pTjet3_lowmass = new TH1D("pTjet3_lowmass", "pT of third leading GenJet", 18, 10., 100.) ;
   mZ->Sumw2(); 
   pTZ->Sumw2();
   nJets->Sumw2();
   pTjet1->Sumw2();
   pTjet2->Sumw2();
   pTjet3->Sumw2();
   
   pTZ_lowmass->Sumw2();
   nJets_lowmass->Sumw2(); 
   pTjet1_lowmass->Sumw2();
   pTjet2_lowmass->Sumw2();
   pTjet3_lowmass->Sumw2();
   return ;
}
 
void LHEPythiaEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;

   edm::Handle< LHEEventProduct > EvtHandle ;
   e.getByLabel( theSrc , EvtHandle ) ;

   float weight = EvtHandle->hepeup().XWGTUP;
   if (whichWeight >= 0) weight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   e.getByLabel(std::string("genParticles"), genParticles);
   TLorentzVector mu[2];
   int iMu = 0;

   for(std::vector<reco::GenParticle>::const_iterator genParticle=genParticles->begin(); genParticle!=genParticles->end(); ++genParticle){
     if( abs(genParticle->pdgId())==13 && genParticle->status()==1 && genParticle->pt() > 10.) {
       mu[iMu].SetPxPyPzE(genParticle->px(),genParticle->py(),
                          genParticle->pz(),genParticle->energy());
       iMu++;
       if (iMu == 2) break;
     }
   }
   TLorentzVector mumu = mu[0] + mu[1];
   
   if (iMu == 2) {
     mZ->Fill(mumu.M(),weight);
     bool lowmass = mumu.M() < 70.;
     if (lowmass) pTZ_lowmass->Fill(mumu.Perp(),weight);
     else pTZ->Fill(mumu.Perp(),weight);
     
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
     
     if (lowmass) {
       nJets_lowmass->Fill(genJetsSorted->size(),weight);
       if (genJetsSorted->size() > 0) pTjet1_lowmass->Fill(genJetsSorted->at(0).pt(),weight);
       if (genJetsSorted->size() > 1) pTjet2_lowmass->Fill(genJetsSorted->at(1).pt(),weight);
       if (genJetsSorted->size() > 2) pTjet3_lowmass->Fill(genJetsSorted->at(2).pt(),weight); 
     } else {
       nJets->Fill(genJetsSorted->size(),weight);
       if (genJetsSorted->size() > 0) pTjet1->Fill(genJetsSorted->at(0).pt(),weight);
       if (genJetsSorted->size() > 1) pTjet2->Fill(genJetsSorted->at(1).pt(),weight);
       if (genJetsSorted->size() > 2) pTjet3->Fill(genJetsSorted->at(2).pt(),weight);
     }
   }
   return ;   
}

void LHEPythiaEventAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(mZ);
  Hlist.Add(pTZ);	   
  Hlist.Add(nJets) ;
  Hlist.Add(pTjet1);
  Hlist.Add(pTjet2); 
  Hlist.Add(pTjet3) ;
  Hlist.Add(pTZ_lowmass);	   
  Hlist.Add(nJets_lowmass) ;
  Hlist.Add(pTjet1_lowmass);
  Hlist.Add(pTjet2_lowmass); 
  Hlist.Add(pTjet3_lowmass) ;
  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n"; 
  cout << "N_events(pT = 0) = " << neventpt0 << "\n"; 
  return ;
}
 
DEFINE_FWK_MODULE(LHEPythiaEventAnalyzer);
