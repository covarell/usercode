// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      HMesonGammaAnalyzer
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
#include "FWCore/Common/interface/TriggerNames.h"

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

#include "Configuration/GenProd/test/HMesonGammaAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
HMesonGammaAnalyzer::HMesonGammaAnalyzer( const ParameterSet& pset )
  : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
    thePhotSrc( pset.getParameter<InputTag>("photonSrc") ), 
    theTrkSrc( pset.getParameter<InputTag>("trackSrc") ),
    thePath( pset.getParameter<string>("HLTriggerName") ),  
    phoChargedIsolationToken_(consumes <edm::ValueMap<float> >
			      (edm::InputTag("photonIDValueMapProducer:phoChargedIsolation")) ),
    phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> >
				    (edm::InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation")) ),
    phoPhotonIsolationToken_(consumes <edm::ValueMap<float> >
			     (edm::InputTag("photonIDValueMapProducer:phoPhotonIsolation")) ),
    rhoToken_(consumes<double> (edm::InputTag("fixedGridRhoFastjetAll")) ),
    trigToken_(consumes<edm::TriggerResults> ( pset.getParameter<InputTag>("HLTriggerResults") )),
    hadronMass( pset.getParameter<double>("hadronMass") ),  
    fOutputFile(0)  
{
   phot_token = mayConsume<edm::View<pat::Photon> >(thePhotSrc);
   trk_token = consumes<std::vector<pat::PackedCandidate> >(theTrkSrc);
}

double HMesonGammaAnalyzer::deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  double deta = eta1 - eta2;
  double dphi = std::fabs(phi1 - phi2);
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return std::sqrt(deta*deta + dphi*dphi);
}

float HMesonGammaAnalyzer::ea(int type, float abseta) { // copyed from Livia Soffi
  if (type == 1) {
    if (abseta < 1) return 0.0385;
    else if (abseta < 1.479) return 0.0468;
    else if (abseta < 2.0) return 0.0435; 
    else if (abseta < 2.2) return 0.0378;
    else if (abseta < 2.3) return 0.0338; 
    else if (abseta < 2.4) return 0.0314; 
    else return 0.0269; 
  }
  if (type == 2) {
    if (abseta < 1) return 0.0636;
    else if (abseta < 1.479) return 0.1103;
    else if (abseta < 2.0) return 0.0759; 
    else if (abseta < 2.2) return 0.0236;
    else if (abseta < 2.3) return 0.0151; 
    else if (abseta < 2.4) return 0.00007; 
    else return 0.0132; 
  }
  if (type == 3) {
    if (abseta < 1) return 0.1240;
    else if (abseta < 1.479) return 0.1093;
    else if (abseta < 2.0) return 0.0631; 
    else if (abseta < 2.2) return 0.0779;
    else if (abseta < 2.3) return 0.0999; 
    else if (abseta < 2.4) return 0.1155; 
    else return 0.1373; 
  }
  return 0.;
}

bool HMesonGammaAnalyzer::trackCuts(pat::PackedCandidate hadron) {
  if (hadron.charge() == 0) return false;
  if (abs(hadron.pdgId()) != 211) return false;
  if (!hadron.trackHighPurity()) return false;
  if (hadron.pt() < 8) return false;
  if (fabs(hadron.eta()) > 2.1) return false;
  return true;
}

void HMesonGammaAnalyzer::beginJob()
{
   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   MTktk = new TH1D("MTktk", "invariant mass of tk-tk",90,0.,4.5) ;
   Minv = new TH1D("Minv", "invariant mass of tk-tk-gamma",80,105.,145.) ;
   Minv_passHLT = new TH1D("Minv_passHLT",  "invariant mass of tk-tk-gamma",80,105.,145.) ;

   nevents = 0; 
   nevPassPho = 0;  
   nevPassIso = 0;
   nevPassMeson = 0;
   return ;
}
 
void HMesonGammaAnalyzer::analyze( const Event& e, const EventSetup& )
{

   nevents++;

   edm::Handle<edm::View<pat::Photon> > photColl;
   e.getByToken( phot_token , photColl ) ;
   edm::Handle< std::vector<pat::PackedCandidate> > trkColl ;
   e.getByToken( trk_token , trkColl ) ;

   edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
   e.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
   edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
   e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
   edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
   e.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
   edm::Handle< double > rhoH;
   e.getByToken(rhoToken_,rhoH);
   rho_ = *rhoH;

   edm::Handle< edm::TriggerResults > trigRes ;
   e.getByToken( trigToken_ , trigRes ) ;
   triggerNames = &( e.triggerNames(*trigRes) );
   unsigned itr = triggerNames->triggerIndex(thePath);

   std::vector<pat::Photon> selPhotons;
   // get just isolated and sorted photons (medium cut-based)
   for (size_t i = 0; i < photColl->size(); ++i){
     const auto pho = photColl->ptrAt(i);
     float pt = pho->pt();
     float abseta = fabs( pho->superCluster()->eta());
     bool isBarrel = (abseta < 1.4442);
     float chIso = (*phoChargedIsolationMap)[pho];
     float nhIso = (*phoNeutralHadronIsolationMap)[pho];
     float phIso = (*phoPhotonIsolationMap)[pho];
     float isoChargedHad = std::max( (float)0.0, chIso - rho_*ea(1,abseta));
     if ((isBarrel && isoChargedHad > 1.416) || (!isBarrel && isoChargedHad > 1.012)) continue;
     float isoNeutralHad = std::max( (float)0.0, nhIso - rho_*ea(2,abseta));
     if ((isBarrel && isoNeutralHad > (2.491 + 0.0126*pt + 0.000026*pt*pt)) || (!isBarrel && isoNeutralHad > (9.131 + 0.0119*pt + 0.000025*pt*pt))) continue;
     float isoPhot = std::max( (float)0.0, phIso - rho_*ea(3,abseta));
     if ((isBarrel && isoPhot > (2.952 + 0.0040*pt)) || (!isBarrel && isoPhot > (4.095 + 0.0040*pt))) continue;
     selPhotons.push_back(*pho);
   }

   if (selPhotons.size() == 0 ) return;
   nevPassPho++;
   std::sort(selPhotons.begin(), selPhotons.end(), pTComparator_);
   pat::Photon thePhoton = selPhotons[0];    // largest pT
   TLorentzVector t1,t2,gamma;

   for (size_t i = 0; i < trkColl->size(); ++i){

     const pat::PackedCandidate trk1 = trkColl->at(i);
     if (!trackCuts(trk1)) continue;

     for (size_t j = i+1; j < trkColl->size(); ++j){

       // find track pairs     
       const pat::PackedCandidate trk2 = trkColl->at(j);
       if (!trackCuts(trk2)) continue;
       if (trk1.charge() * trk2.charge() != -1) continue;
         
       const pat::PackedCandidate trkMaxPt = (trk1.pt() > trk2.pt() ? trk1 : trk2);  
       t1.SetPtEtaPhiM(trk1.pt(),trk1.eta(),trk1.phi(),hadronMass);
       t2.SetPtEtaPhiM(trk2.pt(),trk2.eta(),trk2.phi(),hadronMass);
       TLorentzVector meson = t1 + t2;
       MTktk->Fill(meson.M());

       // rho0 cuts
       if (meson.M() < 0.6 || meson.M() > 0.95) continue;   //rho0 mass
       if (meson.Pt() < 35. || fabs(meson.Rapidity()) > 2.5 ) continue;   //rho0 pT
       nevPassMeson++;

       // rho0 isolation
       float absIso = 0.;
       for (size_t k = 0; k < trkColl->size(); ++k){
	 if (k == i || k == j) continue;
         const pat::PackedCandidate trkN = trkColl->at(k);
	 if (trkN.charge() == 0 || abs(trkN.pdgId()) != 211 || (trkN.fromPV() <= 1 && trkN.dz() > 0.1) || deltaR(trkN.eta(),trkN.phi(),meson.Eta(),meson.Phi()) > 0.3) continue;
	 absIso += trkN.pt();
       }
       if (absIso > 5. && absIso/meson.Pt() > 0.2) continue;
       nevPassIso++;

       // redefine photon eta
       float SCx = thePhoton.superCluster()->x();
       float SCy = thePhoton.superCluster()->y();
       float SCz = thePhoton.superCluster()->z();
       float newTheta = atan(sqrt(SCx*SCx+SCy*SCy)/SCz-trkMaxPt.dz());
       gamma.SetPxPyPzE(thePhoton.energy()*sin(newTheta)*cos(thePhoton.phi()),
			thePhoton.energy()*sin(newTheta)*sin(thePhoton.phi()),				thePhoton.energy()*cos(newTheta),thePhoton.energy() );       
       TLorentzVector higgs = meson + gamma;
       Minv->Fill(higgs.M());
       if ( trigRes->accept(itr) ) Minv_passHLT->Fill(higgs.M());
      
     }
   }
   return ;   
}

void HMesonGammaAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(MTktk);
  Hlist.Add(Minv);
  Hlist.Add(Minv_passHLT);	   
  fOutputFile->cd() ;
  Hlist.Write() ;
  fOutputFile->Close() ;

  std::cout << "Total events = " << nevents << std::endl; 
  std::cout << "Events selected gamma = " << nevPassPho << std::endl;
  std::cout << "Events selected meson = " << nevPassMeson << std::endl;
  std::cout << "Events selected meson isolation = " << nevPassIso << std::endl;
  return ;
}
 
DEFINE_FWK_MODULE(HMesonGammaAnalyzer);
