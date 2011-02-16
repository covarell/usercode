#ifndef HTauTauGenAnalyzer_H
#define HTauTauGenAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"



#include <iostream>
#include <fstream>
 
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
 
// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
 
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObjArray.h"
 
#include "FWCore/Framework/interface/MakerMacros.h"


// forward declarations
class TFile;
class TH1D;
class TH2D;

class HTauTauGenAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit HTauTauGenAnalyzer( const edm::ParameterSet& ) ;
      virtual ~HTauTauGenAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endJob() ;
      int trueVertex(const HepMC::GenEvent* Evt, 
		     const HepMC::GenParticle *aPart);

   private:
   
     //
     std::string fOutputFileName ;
     std::string theSrc ;
     TFile*      fOutputFile ;
     TH1D*       hTauStatus ;	   
     TH1D*       hTauIdDaugs ;
     TH1D*       hPtHiggs ;
     TH1D*       hEtaHiggs ;
     TH1D*       hPtTau ;
     TH1D*       hEtaTau ;
     TH1D*       hPtMu ;
     TH1D*       hEtaMu ;
     TH1D*       hPtEle ;
     TH1D*       hEtaEle ;
     TH1D*       hCosAngEleTau ;
     TH1D*       hCosAngMuTau ; 
     TH1D*       hEtMiss;

     // masses
     TH1D*	 hMtMEtMu ;
     TH1D*	 hMtMEtEle ;
     TH1D*	 hMassMuEle ;
     TH1D*	 hVisibleMass ;
     TH1D*	 hBersaniMass ;
     TH1D*	 hBersaniMassMod ;
     TH1D*	 hSelMtMEtMu ;
     TH1D*	 hSelMtMEtEle ;
     TH1D*	 hSelMassMuEle ;
     TH1D*	 hSelVisibleMass ;
     TH1D*	 hSelBersaniMass ;
     TH1D*	 hSelBersaniMassMod ;
     
     // spins
     TH1D*	 hCosHelAngEle ;     
     TH1D*	 hCosHelAngMu ;	     
     TH1D*	 hCosHelAngTau ;
     TH1D*	 hTransvCosHelAngTau ;
     TH1D*	 hApproxCosHelAngEle ;
     TH1D*	 hApproxCosHelAngMu ;
     TH1D*	 hSelCosHelAngTau ;
     TH1D*	 hSelApproxCosHelAngEle ;
     TH1D*	 hSelApproxCosHelAngMu ;

     // reco electrons/muons/met
     TH1D*       hPtRecoEle ;
     TH1D*       hEtaRecoEle ;
     TH1D*       hMvaRecoEle ;
     TH1D*       hDEtoEtEle ;
     TH1D*       hPtRecoMu ;
     TH1D*       hEtaRecoMu ;
     TH1D*       hDEtoEtMu ;
     TH1D*       hpfMet ;
     TH1D*       hPhipfMet ;
     TH1D*       hpfMetoGenMet ;
     TH2D*       hpfMet_vs_Dr ;

     // reco masses
     TH1D*       hBersaniRecoMass ;
     TH1D*       hSelBersaniRecoMass ;
     TH1D*       hBersaniRecoMassMod ;
     TH1D*       hSelBersaniRecoMassMod ;


     ofstream*   decayed; 
     ofstream*   undecayed; 
     int         nevent, nHiggs;
     
};

#endif
