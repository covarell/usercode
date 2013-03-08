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

#include "GeneratorInterface/ExternalDecays/test/LHEEventAnalyzer.h"

using namespace edm;
using namespace std;
// using namespace HepMC;
 
LHEEventAnalyzer::LHEEventAnalyzer( const ParameterSet& pset )
   : fOutputFileName( pset.getUntrackedParameter<string>("HistOutFile",std::string("HPTRapidity.root")) ),
     theSrc( pset.getUntrackedParameter<string>("theSrc",std::string("source")) ), 
     isVBF( pset.getUntrackedParameter<bool>("isVBF",false) ), 
     fOutputFile(0)
{
}

void LHEEventAnalyzer::beginJob()
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
 
void LHEEventAnalyzer::analyze( const Event& e, const EventSetup& )
{
   nevent++;
   Handle< LHEEventProduct > EvtHandle ;
   
   // find initial HepMCProduct by its label - source
   // OR
   // find HepMCProduct after evtgenlhc by its label - evtgenproducer, that is
   // 
   e.getByLabel( theSrc , EvtHandle ) ;
   
   bool thereIsHiggs = false;

   for (int i = 0; i < EvtHandle->hepeup().NUP; ++i) {
     if (EvtHandle->hepeup().ISTUP[i] != 1) {
       //cout << reader->hepeup.ISTUP[i] << ", " << reader->hepeup.IDUP[i] << endl;
       continue;
     }
     // std::cout << "PDG id LHE = " << EvtHandle->hepeup().IDUP[i] << std::endl;
     if (EvtHandle->hepeup().IDUP[i] == 25) {    // higgs!
       thereIsHiggs = true;
       float mzz = EvtHandle->hepeup().PUP[i][4];
       float pt = sqrt(EvtHandle->hepeup().PUP[i][0]*EvtHandle->hepeup().PUP[i][0] + EvtHandle->hepeup().PUP[i][1]*EvtHandle->hepeup().PUP[i][1]);
       float rapidity = 0.5*log((EvtHandle->hepeup().PUP[i][3] + EvtHandle->hepeup().PUP[i][2])/(EvtHandle->hepeup().PUP[i][3] - EvtHandle->hepeup().PUP[i][2]));
       if (pt < 0.00001) {
	 // std::cout << "pT = " << pt << " mZZ = " << mzz  << std::endl;
	 neventpt0++;
       }
       Pt_sigVBF->Fill(mzz,pt);
       Y_sigVBF->Fill(mzz,rapidity);
     }
   } 

   if (!thereIsHiggs) {   // it's a background event
     
     TLorentzVector z12[2];
     z12[0].SetPxPyPzE(0.,0.,0.,91.18);
     z12[1].SetPxPyPzE(0.,0.,0.,91.18);  // ferm in the rest frame

     int j = 0;
     for (int i = 0; i < EvtHandle->hepeup().NUP; ++i) {
       // std::cout << "PDG id LHE = " << EvtHandle->hepeup().IDUP[i] << std::endl;
       if (EvtHandle->hepeup().IDUP[i] == 23) {    // Z
	 z12[j].SetPxPyPzE(EvtHandle->hepeup().PUP[i][0], 
			   EvtHandle->hepeup().PUP[i][1], 
			   EvtHandle->hepeup().PUP[i][2], 
			   EvtHandle->hepeup().PUP[i][3]);
	 j++; if (j == 2) break;
       
       }
     } 
     
     TLorentzVector tot = z12[0] + z12[1];
     Pt_bkg->Fill(tot.M(),tot.Perp());
     PtOverM_bkg->Fill(tot.M(),tot.Perp()/tot.M());
     Y_bkg->Fill(tot.M(),tot.Rapidity());
   }
   
   return ;   
}

void LHEEventAnalyzer::endJob()
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
 
DEFINE_FWK_MODULE(LHEEventAnalyzer);
