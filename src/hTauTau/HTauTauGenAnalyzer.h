#ifndef HTauTauGenAnalyzer_H
#define HTauTauGenAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

// forward declarations
class TFile;
class TH1D;

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
     
     ofstream*   decayed; 
     ofstream*   undecayed; 
     int         nevent, nHiggs;
     
};

#endif
