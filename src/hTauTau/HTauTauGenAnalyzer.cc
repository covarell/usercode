// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      HTauTauGenAnalyzer
// 
//
// Description: Module to analyze Pythia-EvtGen HepMCproducts
//
//
// Original Author:  Roberto Covarelli
//         Created:  April 26, 2007
//

//#include "GeneratorInterface/ExternalDecays/test/HTauTauGenAnalyzer.h"
#include "UserCode/Covarell/test/HTauTauGenAnalyzer.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"


using namespace edm;
using namespace std;
using namespace HepMC;
 
HTauTauGenAnalyzer::HTauTauGenAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("TestBs.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     fOutputFile(0)
{
}

void HTauTauGenAnalyzer::beginJob()
{
 
   nevent = 0;
   nHiggs = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   
   hTauStatus = new TH1D( "hTauStatus","Status of tau",  5, -0.5, 4.5) ;
   hTauIdDaugs = new TH1D( "hTauIdDaugs","LundIDs of tau daughters",  50, -500., 500.) ;
   hPtHiggs = new TH1D( "hPtHiggs", "Pt Higgs", 50,  0., 170. ) ;
   hPtTau = new TH1D( "hPtTau", "Pt Tau", 50,  0., 170. ) ;
   hPtMu = new TH1D( "hPtMu", "Pt Mu", 50,  0., 100. ) ;
   hPtEle = new TH1D( "hPtEle", "Pt Ele", 50,  0., 100. ) ;
   hEtaHiggs = new TH1D( "hEtaHiggs", "Eta Higgs", 50,  -7.0, 7.0 ) ;
   hEtaTau = new TH1D( "hEtaTau", "Eta Tau", 50,  -7.0, 7.0 ) ;
   hEtaMu = new TH1D( "hEtaMu", "Eta Mu", 50,  -7.0, 7.0 ) ;
   hEtaEle = new TH1D( "hEtaEle", "Eta Ele", 50,  -7.0, 7.0 ) ;
   hCosAngMuTau = new TH1D( "hCosAngMuTau", "cos(#theta_{#tau#mu})", 50,  0.7, 1.0 ) ;
   hCosAngEleTau = new TH1D( "hCosAngEleTau", "cos(#theta_{#tau e})", 50,  0.7, 1.0 ) ;
   hEtMiss = new TH1D( "hEtMiss", "Missing Et", 50,  0., 170. ) ;
   hMtMEtMu = new TH1D( "hMtMEtMu", "Mt EtMiss + Mu", 50,  0., 170. ) ;
   hMtMEtEle = new TH1D( "hMtMEtEle", "Mt EtMiss + Ele", 50,  0., 170. ) ;
   hMassMuEle = new TH1D( "hMassMuEle", "Mass Ele + Mu", 50,  0., 170. ) ;
   hVisibleMass = new TH1D( "hVisibleMass", "Visible Mass", 50,  0., 170. ) ;
   hBersaniMass = new TH1D( "hBersaniMass", "Bersani Mass", 50,  0., 170. ) ;
   hBersaniMassMod = new TH1D( "hBersaniMassMod", "Bersani Mass modified", 50,  0., 170. ) ;
   hSelMtMEtMu = new TH1D( "hSelMtMEtMu", "Mt EtMiss + Mu", 50,  0., 170. ) ;
   hSelMtMEtEle = new TH1D( "hSelMtMEtEle", "Mt EtMiss + Ele", 50,  0., 170. ) ;
   hSelMassMuEle = new TH1D( "hSelMassMuEle", "Mass Ele + Mu", 50,  0., 170. ) ;
   hSelVisibleMass = new TH1D( "hSelVisibleMass", "Visible Mass", 50,  0., 170. ) ;
   hSelBersaniMass = new TH1D( "hSelBersaniMass", "Bersani Mass", 50,  0., 170. ) ;
   hSelBersaniMassMod = new TH1D( "hSelBersaniMassMod", "Bersani Mass modified", 50,  0., 170. ) ;
   hCosHelAngEle = new TH1D( "hCosHelAngEle", "Ele helicity angle", 50,  -1., 1. ) ;     
   hCosHelAngMu = new TH1D( "hCosHelAngMu", "Mu helicity angle", 50,  -1., 1. ) ;     	     
   hCosHelAngTau = new TH1D( "hCosHelAngTau", "Tau helicity angle", 50,  -1., 1. ) ;
   hSelCosHelAngTau = new TH1D( "hSelCosHelAngTau", "Tau helicity angle", 50,  -1., 1. ) ;
   hTransvCosHelAngTau = new TH1D( "hTransvCosHelAngTau", "Transverse tau helicity angle", 50,  -1., 1. ) ; 
   hApproxCosHelAngEle = new TH1D( "hApproxCosHelAngEle", "Tau -> Ele approximated helicity angle", 50,  -1., 1. ) ;     
   hApproxCosHelAngMu = new TH1D( "hApproxCosHelAngMu", "Tau -> Mu approximated helicity angle", 50,  -1., 1. ) ;
   hApproxCosHelAngEle = new TH1D( "hApproxCosHelAngEle", "Tau -> Ele approximated helicity angle", 50,  -1., 1. ) ;     
   hApproxCosHelAngMu = new TH1D( "hApproxCosHelAngMu", "Tau -> Mu approximated helicity angle", 50,  -1., 1. ) ;
   hSelApproxCosHelAngEle = new TH1D( "hSelApproxCosHelAngEle", "Tau -> Ele approximated helicity angle", 50,  -1., 1. ) ;     
   hSelApproxCosHelAngMu = new TH1D( "hSelApproxCosHelAngMu", "Tau -> Mu approximated helicity angle", 50,  -1., 1. ) ;

      //Add RecoQuantities electrons/muons/pfmet
   hEtaRecoEle = new TH1D( "hEtaRecoEle", "Eta Reco Ele", 50,  -4.0, 4.0 ) ;
   hPtRecoEle = new TH1D( "hPtRecoEle", "Pt Reco Ele",  50,  0., 100. ) ;
   hDEtoEtEle = new TH1D( "hDEtoEtEle", "Resol Reco Ele",  50,  -0.5, 0.5 ) ;
   hMvaRecoEle = new TH1D( "hMvaRecoEle", "Resol Reco Ele",  50,  -1., 1. ) ;
   hEtaRecoMu = new TH1D( "hEtaRecoMu", "Eta Reco Mu", 50,  -4.0, 4.0 ) ;
   hPtRecoMu = new TH1D( "hPtRecoMu", "Pt Reco Mu",  50,  0., 100. ) ;
   hDEtoEtMu = new TH1D( "hDEtoEtMu", "Resol Reco Mu",  50,  -0.5, 0.5 ) ;

   // reco pfMET
   hpfMet = new TH1D( "hpfMet", " pf MET ",  50,  0., 100. ) ;
   hPhipfMet = new TH1D( "hPhipfMet", " Phi pf MET ",  50,  -3.2, 3.2 ) ;
   hpfMetoGenMet = new TH1D( "hpfMetoGenMet", "Resol PFMET ",  50,  -1., 1. ) ;
   hpfMet_vs_Dr = new TH2D( "hpfMet_vs_Dr", " PF MET  vs dphi",  50,  0, 3.2, 50, 0., 100. ) ;

   // reco Masses
   hBersaniRecoMass  = new TH1D( "hBersaniRecoMass", "Bersani Reco Mass", 50,  0., 170. ) ;
   hSelBersaniRecoMass  = new TH1D( "hSelBersaniRecoMass", "Bersani Reco Mass", 50,  0., 170. ) ;
   hBersaniRecoMassMod  = new TH1D( "hBersaniRecoMassMod", "Bersani Reco Mass", 50,  0., 170. ) ;
   hSelBersaniRecoMassMod  = new TH1D( "hSelBersaniRecoMassMod", "Bersani Reco Mass", 50,  0., 170. ) ;

   decayed = new ofstream("decayed.txt") ;
   undecayed = new ofstream("undecayed.txt") ;
   return ;
}
 
int HTauTauGenAnalyzer::trueVertex(const GenEvent* Evt, const GenParticle *aPart) {

  // Get rid of fake tau -> tau (gamma) vertices generated by Tauola
  int tauVertId = aPart->end_vertex()->barcode();  
  bool goodTau = false;
  while (!goodTau) {
    goodTau = true;
    const GenVertex* tauvert = Evt->barcode_to_vertex(tauVertId);
    for ( GenVertex::particles_out_const_iterator bp = tauvert->particles_out_const_begin(); bp != tauvert->particles_out_const_end(); ++bp ) {
      if (abs((*bp)->pdg_id()) == 15) {  
	goodTau = false;
	tauVertId = (*bp)->end_vertex()->barcode();
      }
    }
  }
  return tauVertId;

}

void HTauTauGenAnalyzer::analyze( const Event& e, const EventSetup& )
{
      
   Handle< HepMCProduct > EvtHandle ;
   
   // find initial HepMCProduct by its label
   e.getByLabel( theSrc , EvtHandle ) ;
   
   const GenEvent* Evt = EvtHandle->GetEvent() ;
   if (Evt) nevent++;


   // Reco collections
   Handle<reco::PFCandidateCollection> collection;
   InputTag label("particleFlow"); 
   //InputTag label("particleFlow","electrons");  // <- Special electron coll. 
   e.getByLabel(label, collection);
   std::vector<reco::PFCandidate> candidates = (*collection.product());

   Handle<reco::PFMETCollection> pfMETcollection;
   InputTag labelMET("pfMet"); 
   e.getByLabel(labelMET, pfMETcollection);


   for ( GenEvent::particle_const_iterator p = Evt->particles_begin(); p != Evt->particles_end(); ++p ) {

     // look for a Higgs     
     // if ( (*p)->pdg_id() == 25 )  {  // Higgs or Z
     if ( abs((*p)->pdg_id()) > 0 )  {  // any particle  
       GenVertex* endvert = (*p)->end_vertex();
       if (endvert) {

         unsigned int nTau = 0;
	 for ( GenVertex::particles_out_const_iterator ap = endvert->particles_out_const_begin(); ap != endvert->particles_out_const_end(); ++ap ) {
	   if (abs((*ap)->pdg_id()) == 15) nTau++;
	 } 

	 if (nTau != 2) continue;

	 TLorentzVector phiggs((*p)->momentum().px(), 
			       (*p)->momentum().py(),
			       (*p)->momentum().pz(), 
			       (*p)->momentum().e());
	 TVector3 boosterH = - ( phiggs.BoostVector() );

	 TLorentzVector pthiggs((*p)->momentum().px(), 
				(*p)->momentum().py(), 0., 
				sqrt(8315.61 + pow((*p)->momentum().px(),2) 
				     + pow((*p)->momentum().py(),2)));
	 TVector3 boosterHt = - ( pthiggs.BoostVector() );

	 hPtHiggs->Fill((*p)->momentum().perp());
	 hEtaHiggs->Fill((*p)->momentum().pseudoRapidity());

	 // loop 1 --> look for H -> tau tau -> e mu 
	 int theCodeProduct = 1;
	 TLorentzVector pmiss;
	 TLorentzVector pmu;
	 TLorentzVector pele;
	 float theTrueHelAngle = 0.;

	 TLorentzVector precomiss;
	 TLorentzVector precomu;
	 TLorentzVector precoele;
	 bool isThereRecoEle = false;
	 bool isThereRecoMu = false;

	 for ( GenVertex::particles_out_const_iterator ap = endvert->particles_out_const_begin(); ap != endvert->particles_out_const_end(); ++ap ) {
	   if (abs((*ap)->pdg_id()) == 15) { 

	     hPtTau->Fill((*ap)->momentum().perp());
	     hTauStatus->Fill((*ap)->status());
	     hEtaTau->Fill((*ap)->momentum().pseudoRapidity());
	     *decayed << (*ap)->pdg_id() << " --> ";

	     TLorentzVector pta((*ap)->momentum().px(), 
				(*ap)->momentum().py(),
				(*ap)->momentum().pz(), 
				(*ap)->momentum().e());
	     TLorentzVector ptta((*ap)->momentum().px(), 
				 (*ap)->momentum().py(), 0., 
				 sqrt(3.1577 + pow((*p)->momentum().px(),2) 
				      + pow((*p)->momentum().py(),2)));
	     pta.Boost( boosterH );
	     theTrueHelAngle = cos( pta.Vect().Angle(phiggs.Vect()));
	     hCosHelAngTau->Fill( theTrueHelAngle );
	     ptta.Boost( boosterHt );
	     hTransvCosHelAngTau->Fill( cos( ptta.Vect().Angle(pthiggs.Vect())) );
	     
	     int tauVertId = trueVertex(Evt, *ap);  
	    
             const GenVertex* tauvert2 = Evt->barcode_to_vertex(tauVertId);
	     for ( GenVertex::particles_out_const_iterator cp = tauvert2->particles_out_const_begin(); cp != tauvert2->particles_out_const_end(); ++cp ) {
	       
	       hTauIdDaugs->Fill((*cp)->pdg_id());
	       *decayed << (*cp)->pdg_id() << " ";
	       if (abs((*cp)->pdg_id()) == 11 || abs((*cp)->pdg_id()) == 13) 
		 theCodeProduct *= (*cp)->pdg_id();
	     }
	     *decayed << "\n";
	   }
	 }

	 if (theCodeProduct != -143) continue; // only e-mu opposite sign 
	                                       // -(13*11) = -143
	 
	 nHiggs++;
       
	 // loop 2 --> fill e mu quantities 
	 for ( GenVertex::particles_out_const_iterator ap = endvert->particles_out_const_begin(); ap != endvert->particles_out_const_end(); ++ap ) {
	   if (abs((*ap)->pdg_id()) == 15) {
	      
	     int tauVertId = trueVertex(Evt, *ap);

             const GenVertex* tauvert2 = Evt->barcode_to_vertex(tauVertId);
	     for ( GenVertex::particles_out_const_iterator cp = tauvert2->particles_out_const_begin(); cp != tauvert2->particles_out_const_end(); ++cp ) {

	       TLorentzVector pta((*ap)->momentum().px(), 
				  (*ap)->momentum().py(),
				  (*ap)->momentum().pz(), 
				  (*ap)->momentum().e());
	       TVector3 boosterTa = - ( pta.BoostVector() );
	       TLorentzVector pl((*cp)->momentum().px(), 
				 (*cp)->momentum().py(),
				 (*cp)->momentum().pz(), 
				 (*cp)->momentum().e());
	       if (abs((*cp)->pdg_id()) == 11) {  // e
		 pele = pl;
		 // lab frame
		 hPtEle->Fill((*cp)->momentum().perp());
		 hEtaEle->Fill((*cp)->momentum().pseudoRapidity());     
		 hCosAngEleTau->Fill( cos( pl.Vect().Angle(pta.Vect())));
		 // own frames
		 pl.Boost( boosterTa );
		 hCosHelAngEle->Fill( cos( pl.Vect().Angle(pta.Vect())) );

		 	 // find the electrons matched to generated.
		 std::vector<reco::PFCandidate>::iterator it;
		 for ( it = candidates.begin(); it != candidates.end(); ++it)   {
		   
		   reco::PFCandidate::ParticleType type = (*it).particleId();
		   if ( type == reco::PFCandidate::e) {
		     GsfTrackRef refGsf = (*it).gsfTrackRef();
		     unsigned int nLostHits = 
		       refGsf->trackerExpectedHitsInner().numberOfLostHits();
		     // reject conversions
		     

		     if(nLostHits == 0) {
		       // compute electron isolation
		       float pxcand = (*it).px();
		       float pycand = (*it).py();
		       float pzcand = (*it).pz();
		       float enecand =  (*it).energy();
		       
		       TLorentzVector cand_ele;
		       cand_ele.SetPxPyPzE(pxcand,pycand,pzcand,enecand);
		       
		       float deta = fabs(cand_ele.Eta() - pele.Eta());
		       float dphi =  fabs(cand_ele.Phi() - pele.Phi());
		       if (dphi>TMath::Pi()) dphi-= TMath::TwoPi();
		       float dR = sqrt(deta*deta +
				       dphi*dphi);
		       if(dR < 0.1) {
			 hPtRecoEle->Fill(cand_ele.Perp());   
 			 hEtaRecoEle->Fill(cand_ele.Eta());   
			 float DEtoEt =  (cand_ele.Perp() - pele.Perp()) / pele.Perp();
			 hDEtoEtEle->Fill(DEtoEt);
 			 hMvaRecoEle->Fill((*it).mva_e_pi());
			 // add isolation and add GsfElectrons
			 precoele = cand_ele;
			 isThereRecoEle = true;
		       }
		       
		     }
		   }
		 } // end loop on candidates
		 
		 


	       } else if (abs((*cp)->pdg_id()) == 13) {  // mu
		 pmu = pl;
		 // lab frame
		 hPtMu->Fill((*cp)->momentum().perp());
		 hEtaMu->Fill((*cp)->momentum().pseudoRapidity());
		 hCosAngMuTau->Fill( cos( pl.Vect().Angle(pta.Vect()))) ;
		 // own frames
		 pl.Boost( boosterTa );
		 hCosHelAngMu->Fill( cos( pl.Vect().Angle(pta.Vect())) );

		 std::vector<reco::PFCandidate>::iterator it;
		 for ( it = candidates.begin(); it != candidates.end(); ++it)   {
		   
		   reco::PFCandidate::ParticleType type = (*it).particleId();
		   if ( type == reco::PFCandidate::mu) {
		     float pxcand = (*it).px();
		     float pycand = (*it).py();
		     float pzcand = (*it).pz();
		     float enecand =  (*it).energy();
		     
		     TLorentzVector cand_mu;
		     cand_mu.SetPxPyPzE(pxcand,pycand,pzcand,enecand);
		     
		     float deta = fabs(cand_mu.Eta() - pmu.Eta());
		     float dphi =  fabs(cand_mu.Phi() - pmu.Phi());
		     if (dphi>TMath::Pi()) dphi-= TMath::TwoPi();
		     float dR = sqrt(deta*deta +
				     dphi*dphi);
		     if(dR < 0.1) {
		       hPtRecoMu->Fill(cand_mu.Perp());   
		       hEtaRecoMu->Fill(cand_mu.Eta());   
		       float DEtoEt =  (cand_mu.Perp() - pmu.Perp()) / pmu.Perp();
		       hDEtoEtMu->Fill(DEtoEt);
		       //  hMvaRecoMu->Fill((*it).mva_e_pi());
		       // add isolation and add GsfMuctrons
		       precomu = cand_mu;
		       isThereRecoMu = true;
		     }
		     
		   }
		 } // end loop on candidates


	       } else {   // neutrinos
		 pmiss += pl;
	       }
	     }
	   }
	 }
	 
	 hEtMiss->Fill(pmiss.Perp());

         // Calculate h -> tau tau observables
	 TLorentzVector pteleEst(pele.X(),pele.Y(),
			      0.,sqrt(pow(pele.X(),2) + pow(pele.Y(),2)));
	 TLorentzVector ptmuEst(pmu.X(),pmu.Y(),
			      0.,sqrt(pow(pmu.X(),2) + pow(pmu.Y(),2)));
	 TLorentzVector ptmissEst(pmiss.X(),pmiss.Y(),
			      0.,sqrt(pow(pmiss.X(),2) + pow(pmiss.Y(),2)));
	 TLorentzVector emetEst = pteleEst + ptmissEst;
	 TLorentzVector mumetEst = ptmuEst + ptmissEst;
	 TLorentzVector emu = pele + pmu;
	 TLorentzVector emumissEst = emu + ptmissEst;
         TLorentzVector pthiggsEst = pteleEst + ptmuEst + ptmissEst;

	 // Bersani calculation
	 float xtaue = (pele.X()*pmu.Y() - pele.Y()*pmu.X())/(pthiggs.X()*pmu.Y() - pthiggs.Y()*pmu.X());
	 float xtaumu = (pele.Y()*pmu.X() - pele.X()*pmu.Y())/(pthiggs.X()*pele.Y() - pthiggs.Y()*pele.X());
	 // cout << xtau1 << " " << xtau2 << endl;

         if (xtaue*xtaumu < 0. || fabs(xtaue) > 2. || fabs(xtaumu) > 2.) {
	   cout << endl << "Estimated tau momentum fraction is negative or much greater than 1." << endl;
	   cout << "x_taue = " << xtaue << " x_taumu = " << xtaumu << endl;
	   cout << "Skipping ... " << endl;
	   continue;
	 }

	 TLorentzVector ptaueEst = pele*(1./xtaue);
         TLorentzVector ptaumuEst = pmu*(1./xtaumu);
	 // force tau mass - probably not useful
	 ptaueEst.SetE(sqrt(3.1577 + pow(ptaueEst.X(),2) + 
			    pow(ptaueEst.Y(),2) + pow(ptaueEst.Z(),2)) );
	 ptaumuEst.SetE(sqrt(3.1577 + pow(ptaumuEst.X(),2) + 
			     pow(ptaumuEst.Y(),2) + pow(ptaumuEst.Z(),2)) );
	 TLorentzVector phiggsEst = ptaueEst + ptaumuEst;
         TVector3 boosterHEst = - ( phiggsEst.BoostVector() );
	 
	 hMtMEtEle->Fill(emetEst.Perp());
	 hMtMEtMu->Fill(mumetEst.Perp());
	 hMassMuEle->Fill(emu.M());
	 hVisibleMass->Fill(emumissEst.M());
	 hBersaniMass->Fill(emu.M()/sqrt(xtaue*xtaumu));
	 hBersaniMassMod->Fill(phiggsEst.M());    
	 ptaueEst.Boost( boosterHEst );
	 ptaumuEst.Boost( boosterHEst );
	 hApproxCosHelAngEle->Fill( cos( ptaueEst.Vect().Angle(phiggsEst.Vect())) );
	 hApproxCosHelAngMu->Fill( cos( ptaumuEst.Vect().Angle(phiggsEst.Vect())) );
	 
	 // Reasonable offline cuts
	 if (pmu.Perp() > 10 && pele.Perp() > 15) {
	   hSelCosHelAngTau->Fill( theTrueHelAngle );
	   hSelMtMEtEle->Fill(emetEst.Perp());
	   hSelMtMEtMu->Fill(mumetEst.Perp());
	   hSelMassMuEle->Fill(emu.M());
	   hSelVisibleMass->Fill(emumissEst.M());
	   hSelBersaniMass->Fill(emu.M()/sqrt(xtaue*xtaumu));
	   hSelBersaniMassMod->Fill(phiggsEst.M());
	   hSelApproxCosHelAngEle->Fill( cos( ptaueEst.Vect().Angle(pthiggsEst.Vect())) );
	   hSelApproxCosHelAngMu->Fill( cos( ptaumuEst.Vect().Angle(pthiggsEst.Vect())) );

	 }



	 // Do the same with reconstructed variables


	 if(isThereRecoEle && isThereRecoMu) {
	   float pfMet = pfMETcollection->begin()->et();
	   float pfMet_px = pfMETcollection->begin()->px();
	   float pfMet_py = pfMETcollection->begin()->py();
	   float pfMet_phi = pfMETcollection->begin()->phi();
	   hpfMet->Fill(pfMet);
	   hPhipfMet->Fill(pfMet_phi);
	   
	   //resolution of the pfMET
	   hpfMetoGenMet->Fill(pfMet);
	   float lepDphi = fabs(precoele.Phi() - precomu.Phi());
	   float lepDeta = fabs(precoele.Eta() - precomu.Eta());
	   if (lepDphi>TMath::Pi()) lepDphi-= TMath::TwoPi();
	   float lepdR = sqrt(lepDeta*lepDeta +
			   lepDphi*lepDphi);

	   hpfMet_vs_Dr->Fill(lepdR,pfMet);



	   if(pfMet > 20) {
	     
	     cout << " Gen MET px " << pmiss.X() 
		  << " Gen MET py " << pmiss.Y() 
		  << " Gen MET    " << pmiss.Perp() << endl;
	     
	     cout << " Reco MET px " <<  pfMet_px
		  << " Reco MET py " <<  pfMet_py
		  << " Reco MET    " <<  pfMet   << endl;
	     
	     
	     //	 float ptrecomiss = pfMet;
	     TLorentzVector ptrecoele(precoele.X(),precoele.Y(),
				      0.,sqrt(pow(precoele.X(),2) + pow(precoele.Y(),2)));
	     TLorentzVector ptrecomu(precomu.X(),precomu.Y(),
				     0.,sqrt(pow(precomu.X(),2) + pow(precomu.Y(),2)));
	     TLorentzVector ptrecomiss(pfMet_px,pfMet_py,
				       0.,sqrt(pow(pfMet_px,2) + pow(pfMet_py,2)));
	     
	     TLorentzVector recoemu = precoele + precomu;
	     TLorentzVector ptrecohiggs = ptrecoele + ptrecomu + ptrecomiss;
	     float xrecotaue = 
	       (precoele.X()*precomu.Y() - precoele.Y()*precomu.X())/(ptrecohiggs.X()*precomu.Y() - ptrecohiggs.Y()*precomu.X());
	     float xrecotaumu = 
	       (precoele.Y()*precomu.X() - precoele.X()*precomu.Y())/(ptrecohiggs.X()*precoele.Y() - ptrecohiggs.Y()*precoele.X());


	     if (xrecotaue*xrecotaumu < 0. || fabs(xrecotaue) > 2. || fabs(xrecotaumu) > 2.) {
	       cout << endl << "Estimated RECO tau momentum fraction is negative or much greater than 1." << endl;
	       cout << "x_taue = " << xrecotaue << " x_taumu = " << xrecotaumu << endl;
	       cout << "Skipping RECO... " << endl;
	       continue;
	     }
	     
	     TLorentzVector precotaueEst = precoele*(1./xrecotaue);
	     TLorentzVector precotaumuEst = precomu*(1./xrecotaumu);
	     precotaueEst.SetE(sqrt(3.1577 + pow(precotaueEst.X(),2) + 
				    pow(precotaueEst.Y(),2) + pow(precotaueEst.Z(),2)) );
	     precotaumuEst.SetE(sqrt(3.1577 + pow(precotaumuEst.X(),2) + 
				     pow(precotaumuEst.Y(),2) + pow(precotaumuEst.Z(),2)) );
	     
	     TLorentzVector precohiggsEst = precotaueEst + precotaumuEst;
	     TVector3 boosterRecoHEst = - ( precohiggsEst.BoostVector() );
	     


	     hBersaniRecoMass->Fill(recoemu.M()/sqrt(xrecotaue*xrecotaumu));
	     hBersaniRecoMassMod->Fill(precohiggsEst.M());

	     if (precomu.Perp() > 10 && precoele.Perp() > 15) {
	       hSelBersaniRecoMass->Fill(recoemu.M()/sqrt(xrecotaue*xrecotaumu));
	       hSelBersaniRecoMassMod->Fill(precohiggsEst.M());
	     }

	     


	   }
	 }
       }
     }
   }

   return ;   
}

void HTauTauGenAnalyzer::endJob()
{
  TObjArray Hlist(0);
  Hlist.Add(hTauStatus) ;	   
  Hlist.Add(hTauIdDaugs) ;
  Hlist.Add(hPtHiggs) ;
  Hlist.Add(hEtaHiggs) ;
  Hlist.Add(hPtTau) ;
  Hlist.Add(hEtaTau) ;
  Hlist.Add(hPtMu) ;
  Hlist.Add(hEtaMu) ;
  Hlist.Add(hPtEle) ;
  Hlist.Add(hEtaEle) ;
  Hlist.Add(hCosAngMuTau) ;
  Hlist.Add(hCosAngEleTau) ;
  Hlist.Add(hEtMiss);
  Hlist.Add(hMtMEtMu); 
  Hlist.Add(hMtMEtEle);
  Hlist.Add(hMassMuEle);
  Hlist.Add(hVisibleMass);
  Hlist.Add(hBersaniMass);
  Hlist.Add(hBersaniMassMod);
  Hlist.Add(hSelMtMEtMu); 
  Hlist.Add(hSelMtMEtEle);
  Hlist.Add(hSelMassMuEle);
  Hlist.Add(hSelVisibleMass);
  Hlist.Add(hSelBersaniMass);
  Hlist.Add(hSelBersaniMassMod);
  Hlist.Add(hCosHelAngEle) ;       
  Hlist.Add(hCosHelAngMu) ;	     
  Hlist.Add(hCosHelAngTau) ;  
  Hlist.Add(hSelCosHelAngTau) ;  
  Hlist.Add(hTransvCosHelAngTau) ;
  Hlist.Add(hApproxCosHelAngEle) ;
  Hlist.Add(hApproxCosHelAngMu) ;
  Hlist.Add(hSelApproxCosHelAngEle) ;
  Hlist.Add(hSelApproxCosHelAngMu) ;
  Hlist.Add(hPtRecoEle) ;
  Hlist.Add(hEtaRecoEle) ;
  Hlist.Add(hMvaRecoEle) ;
  Hlist.Add(hDEtoEtEle) ;
  Hlist.Add(hPtRecoMu) ;
  Hlist.Add(hEtaRecoMu) ;
  Hlist.Add(hDEtoEtMu) ;
  Hlist.Add(hpfMet) ;
  Hlist.Add(hPhipfMet) ;
  Hlist.Add(hpfMetoGenMet) ;
  Hlist.Add(hBersaniRecoMass) ;
  Hlist.Add(hSelBersaniRecoMass) ;
  Hlist.Add(hpfMet_vs_Dr) ;
  Hlist.Add(hBersaniRecoMassMod) ;
  Hlist.Add(hSelBersaniRecoMassMod) ;

  //  Hlist.Add() ;

  Hlist.Write() ;
  fOutputFile->Close() ;
  cout << "N_events = " << nevent << "\n";
  cout << "N(H -> tau tau -> e mu) = " << nHiggs << "\n"; 
  return ;
}
 
DEFINE_FWK_MODULE(HTauTauGenAnalyzer);
