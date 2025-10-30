// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      HPhiPhiAnalyzer
// 
//
// Description: Module to analyze H -> rho/phi gamma 
//
//
// Original Author:  Roberto Covarelli
//         Created:  May 29, 2018
//

#include <iostream>
#include <fstream>
 
// #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObjArray.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Hrhogamma/Hrhogamma/test/HPhiPhiAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
HPhiPhiAnalyzer::HPhiPhiAnalyzer( const ParameterSet& pset )
  : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
    theTrkSrc( pset.getParameter<InputTag>("trackSrc") ),
    theTrkSrc2( pset.getParameter<InputTag>("trackSrc2") ),
    thePath( pset.getParameter<string>("HLTriggerName") ),  
    trigToken_(consumes<edm::TriggerResults> ( pset.getParameter<InputTag>("HLTriggerResults") )),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(pset.getParameter<InputTag>("HLTriggerObjects"))),
    hadronMass( pset.getParameter<double>("hadronMass") ),  
    fOutputFile(0)  
{
  // all initializations
  trk_token = consumes<std::vector<pat::PackedCandidate> >(theTrkSrc);
  trk_token2 = consumes<std::vector<pat::PackedCandidate> >(theTrkSrc2);
}

double HPhiPhiAnalyzer::deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  double deta = eta1 - eta2;
  double dphi = std::fabs(phi1 - phi2);
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return std::sqrt(deta*deta + dphi*dphi);
}

bool HPhiPhiAnalyzer::trackCuts(pat::PackedCandidate hadron) {
  if (hadron.charge() == 0) return false;
  if (abs(hadron.pdgId()) != 211) return false;
  if (!hadron.trackHighPurity()) return false;
  if (hadron.pt() < 8.) return false;
  if (fabs(hadron.eta()) > 2.4) return false;
  return true;
}

void HPhiPhiAnalyzer::beginJob()
{
   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // invariant masses
   MTktk1 = new TH1D("MTktk1", "invariant mass of tk-tk",90,0.,4.5) ;
   MTktk2 = new TH1D("MTktk2", "invariant mass of tk-tk",90,0.,4.5) ;
   Minv = new TH1D("Minv", "invariant mass of 4tk",80,105.,145.) ;
   Minv_passHLT = new TH1D("Minv_passHLT",  "invariant mass of 4tk",80,105.,145.) ;

   // pt-eta
   pTmeson = new TH1D("pTmeson", "pT of meson",20,15.,200.) ;
   etameson = new TH1D("etameson", "pT of gamma",20,-2.6,2.6) ;
       
   // trigger efficiencies vs. pT-eta
   effpTmeson = new TEfficiency("effpTmeson","my efficiency;pT;#epsilon",20,15.,200.);
   effetameson = new TEfficiency("effetameson","my efficiency;eta;#epsilon",20,-2.6,2.6);

   // counters
   nevents = 0; 
   neventsTauMatch = 0;
   neventsTauPass = 0;
   return ;
}
 
void HPhiPhiAnalyzer::analyze( const Event& e, const EventSetup& )
{

   nevents++;

   // get all objects
   edm::Handle< std::vector<pat::PackedCandidate> > trkColl ;
   e.getByToken( trk_token , trkColl ) ;
   edm::Handle< std::vector<pat::PackedCandidate> > trkColl2 ;
   e.getByToken( trk_token2 , trkColl2 ) ;

   edm::Handle< edm::TriggerResults > trigRes ;
   e.getByToken( trigToken_ , trigRes ) ;
   edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
   e.getByToken(triggerObjects_, triggerObjects);
   triggerNames = &( e.triggerNames(*trigRes) );
   unsigned itr = triggerNames->triggerIndex(thePath);

   // Merge tracks and "lost tracks" (???)
   std::vector<pat::PackedCandidate> allTracks;
   for (size_t i = 0; i < trkColl->size(); ++i){
     const pat::PackedCandidate trk = trkColl->at(i);
     allTracks.push_back(trk);
   }
   for (size_t j = 0; j < trkColl2->size(); ++j){
     const pat::PackedCandidate trk2 = trkColl2->at(j);
     allTracks.push_back(trk2);
   }

   TLorentzVector t1,t2;
   TLorentzVector theMeson[2];
   float theMesonPt = 0.;
   float theDz = 0.;

   // find track pairs 
   for (size_t i = 0; i < allTracks.size(); ++i){

     const pat::PackedCandidate trk1 = allTracks.at(i);
     if (!trackCuts(trk1)) continue;

     for (size_t j = i+1; j < allTracks.size(); ++j){
   
       const pat::PackedCandidate trk2 = allTracks.at(j);
       if (!trackCuts(trk2)) continue;
       if (trk1.charge() * trk2.charge() != -1) continue;
       if (fabs(trk1.pt() - trk2.pt()) < 0.01 && fabs(trk1.phi() - trk2.phi()) < 0.01 && fabs(trk1.eta() - trk2.eta()) < 0.01) continue;  // avoid duplication
       if deltaR(trk1.eta(),trk2.phi(),trk1.Eta(),trk2.Phi()) > 0.5) continue;
         
       const pat::PackedCandidate trkMaxPt = (trk1.pt() > trk2.pt() ? trk1 : trk2);  
       t1.SetPtEtaPhiM(trk1.pt(),trk1.eta(),trk1.phi(),hadronMass);
       t2.SetPtEtaPhiM(trk2.pt(),trk2.eta(),trk2.phi(),hadronMass);
       TLorentzVector meson = t1 + t2;
       MTktk->Fill(meson.M());

       // rho0 cuts
       if (meson.M() < 0.5 || meson.M() > 1.0) continue;   //rho0 mass
       if (meson.Pt() < 35. || fabs(meson.Rapidity()) > 2.1 ) continue;   //rho0 pT
     
       // rho0 isolation
       float absIso = 0.;
       for (size_t k = 0; k < allTracks.size(); ++k){
	 if (k == i || k == j) continue;
         const pat::PackedCandidate trkN = allTracks.at(k);
	 if (trkN.charge() == 0 || trkN.pt() < 0.5 ||
	     abs(trkN.pdgId()) != 211 || 
	     ( trkN.fromPV() <= 1 && trkN.dz() > 0.1 ) || 
	     deltaR(trkN.eta(),trkN.phi(),meson.Eta(),meson.Phi()) > 0.5) continue;
	 absIso += trkN.pt();
       }
       if (absIso/meson.Pt() > 0.2) continue;
       
       if (meson.Pt() > theMesonPt) {
	 theMesonPt = meson.Pt();
	 theMeson = meson;             // largest pT
         theDz = trkMaxPt.dz();
       }
     }
   }

   if (theMesonPt < 1.) return;
   neventsTauPass++;
   pTmeson->Fill(theMeson.Pt());
   etameson->Fill(theMeson.Eta());

   // trigger matches
   bool thisTauMatch = false;
   for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      bool lastFilterOK = false;
      obj.unpackPathNames(*triggerNames);
      obj.unpackFilterLabels(e,*trigRes);
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) if (obj.filterLabels()[h] == "hltDoublePFJets30PNetTauhTagMediumWP") lastFilterOK = true;
      if (!lastFilterOK) continue;
      if (obj.collection() == "hltPFJetForBtag::HLT" && deltaR(theMeson.Eta(),theMeson.Phi(),obj.eta(),obj.phi()) < 0.3) {
	neventsTauMatch++;  thisTauMatch = true;  break;
      }
   }	
   effpTmeson->Fill(thisTauMatch,theMeson.Pt());
   effetameson->Fill(thisTauMatch,theMeson.Eta());

   TLorentzVector higgs = theMeson + gamma;
   Minv->Fill(higgs.M());
   if ( trigRes->accept(itr) ) Minv_passHLT->Fill(higgs.M());

   return ;   
}

void HPhiPhiAnalyzer::endJob()
{

  TObjArray Hlist(0);
  Hlist.Add(MTktk1);
  Hlist.Add(MTktk2);
  Hlist.Add(Minv);
  Hlist.Add(Minv_passHLT);
  Hlist.Add(pTmeson);
  Hlist.Add(etameson);
  Hlist.Add(effpTmeson);
  Hlist.Add(effetameson);

  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;

  std::cout << "Total events = " << nevents << std::endl; 
  std::cout << "Events selected gamma+tau = " << neventsTauPass << std::endl;
  std::cout << "Events matched gamma+tau = " << neventsTauMatch << std::endl;
  return ;
}
 
DEFINE_FWK_MODULE(HPhiPhiAnalyzer);
