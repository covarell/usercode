// -*- C++ -*-
//
// Package:    GeneratorInterface
// Class:      BsAnalyzer
// 
//
// Description: Module to analyze Pythia-EvtGen HepMCproducts
//
//
// Original Author:  Roberto Covarelli
//         Created:  April 26, 2007
//

#include <iostream>
#include "fstream.h"
 
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
 
// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"

#include "GeneratorInterface/EvtGenInterface/test/BsAnalyzer.h"

using namespace edm;
using namespace std;

 
BsAnalyzer::BsAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("TestBs.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     fOutputFile(0)
{
}

void BsAnalyzer::beginJob( const EventSetup& )
{
 
   nevent = 0;
   nbs = 0;

   fOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   // fHist2muMass  = new TH1D(  "Hist2muMass", "2-mu inv. mass", 100,  60., 120. ) ;  
   hGeneralId = new TH1D( "hGeneralId","LundIDs of all particles",  100, -1000., 1000.) ;
   hPtbs = new TH1D( "hPtbs", "Pt Bs", 100,  0., 50. ) ;
   hPbs  = new TH1D( "hPbs",  "P Bs",  100,  0., 200. ) ;
   hPhibs = new TH1D( "hPhibs","Phi Bs",  100,  -3.14, 3.14) ;
   hEtabs = new TH1D( "hEtabs","Eta Bs",  100,  -15.0, 15.0) ;
   hPtmu = new TH1D( "hPtmu", "Pt Mu", 100,  0., 50. ) ;
   hPmu  = new TH1D( "hPmu",  "P Mu",  100,  0., 200. ) ;
   hPhimu = new TH1D( "hPhimu","Phi Mu",  100,  -3.14, 3.14) ;
   hEtamu = new TH1D( "hEtamu","Eta Mu",  100,  -15.0, 15.0) ;
   hmumuMass = new TH1D( "hmumuMass","#mu^{+}#mu^{-} invariant mass",  100, -1.0, 10.0) ;
   hIdBsDaugs = new TH1D( "hIdBsDaugs","LundIDs of the Bs's daughters",  100, -1000., 1000.) ;
   hIdPhiDaugs = new TH1D( "hIdPhiDaugs","LundIDs of the phi's daughters",  100, -500., 500.) ;
   hIdBDaugs = new TH1D( "hIdBDaugs","LundIDs of the B's daughters",  100, -1000., 1000.) ;

   decayed = new ofstream("decayed.txt") ;
   undecayed = new ofstream("undecayed.txt") ;
   return ;
}
 
void BsAnalyzer::analyze( const Event& e, const EventSetup& )
{
      
   Handle< HepMCProduct > EvtHandle ;
   
   // find initial HepMCProduct by its label - source
   // OR
   // find HepMCProduct after evtgenlhc by its label - evtgenproducer, that is
   // 
   e.getByLabel( theSrc , EvtHandle ) ;
   
   const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
   if (Evt) nevent++;

   for ( HepMC::GenEvent::particle_const_iterator p = Evt->particles_begin(); p != Evt->particles_end(); ++p ) {
     HepMC::GenVertex* genvert = (*p)->end_vertex(); 
     hGeneralId->Fill((*p)->pdg_id()); 
     if ( abs((*p)->pdg_id()/100) == 5 || abs((*p)->pdg_id()/100) == 4 || abs((*p)->pdg_id()/100) == 3) {
       if (!genvert) {
	 *undecayed << (*p)->pdg_id() << endl;
       } else {
	 *decayed << (*p)->pdg_id() << endl;
       }
     }
     // if ( abs((*p)->pdg_id()) == 211 || abs((*p)->pdg_id()) == 213 ) (*p)->print(); // pi vs rho
     if ( (*p)->pdg_id() == 531 ) {  // B_s 
       nbs++;
       hPtbs->Fill((*p)->momentum().perp());
       hPbs->Fill( sqrt ( pow((*p)->momentum().px(),2)+pow((*p)->momentum().py(),2)+
			  pow((*p)->momentum().pz(),2) )) ;
       hPhibs->Fill((*p)->momentum().phi());
       hEtabs->Fill((*p)->momentum().pseudoRapidity());
       if (genvert) {
         for ( HepMC::GenVertex::particles_out_const_iterator ap = genvert->particles_out_const_begin(); ap != genvert->particles_out_const_end(); ++ap ) {
	   hIdBsDaugs->Fill((*ap)->pdg_id());
	 }
       }
     }
     if ( abs((*p)->pdg_id()) == 511 ) {  // B0 
       if (genvert) {
         for ( HepMC::GenVertex::particles_out_const_iterator bp = genvert->particles_out_const_begin(); bp != genvert->particles_out_const_end(); ++bp ) {
	   hIdBDaugs->Fill((*bp)->pdg_id());
	 }
       }
     }
     if ( (*p)->pdg_id() == 333 ) {  // phi 
       if (genvert) {
         for ( HepMC::GenVertex::particles_out_const_iterator cp = genvert->particles_out_const_begin(); cp != genvert->particles_out_const_end(); ++cp ) {
	   hIdPhiDaugs->Fill((*cp)->pdg_id());
	 }
       }
     }
     if ( (*p)->pdg_id() == 13 ) { // mu+
       hPtmu->Fill((*p)->momentum().perp());
       hPmu->Fill( sqrt ( pow((*p)->momentum().px(),2)+pow((*p)->momentum().py(),2)+
			  pow((*p)->momentum().pz(),2) )) ;
       hPhimu->Fill((*p)->momentum().phi());
       hEtamu->Fill((*p)->momentum().pseudoRapidity());
       for ( HepMC::GenEvent::particle_const_iterator p2 = Evt->particles_begin(); p2 != Evt->particles_end(); ++p2 ) {
	 if ( (*p2)->pdg_id() == -13 ) { // mu-
	   HepMC::FourVector pdiff;
           cout << "MU+ " << (*p)->momentum().px() << " " << (*p)->momentum().py() << " " << (*p)->momentum().pz() << " " << (*p)->momentum().e() << endl;
           cout << "MU- " << (*p2)->momentum().px() << " " << (*p2)->momentum().py() << " " << (*p2)->momentum().pz() << " " << (*p2)->momentum().e() << endl;
           pdiff.set((*p)->momentum().px() + (*p2)->momentum().px(),
		     (*p)->momentum().py() + (*p2)->momentum().py(),
                     (*p)->momentum().pz() + (*p2)->momentum().pz(),
		     (*p)->momentum().e() + (*p2)->momentum().e());
           hmumuMass->Fill(pdiff.m());
	 }
       }
     }
   }

   return ;
   
}

void BsAnalyzer::endJob()
{
       
   fOutputFile->Write() ;
   fOutputFile->Close() ;
   cout << "N_events = " << nevent << "\n";
   cout << "N_Bs = " << nbs << "\n"; 
   return ;
}
 
DEFINE_FWK_MODULE(BsAnalyzer);
