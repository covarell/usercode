#ifndef BsAnalyzer_H
#define BsAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

// forward declarations
class TFile;
class TH1D;

class BsAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit BsAnalyzer( const edm::ParameterSet& ) ;
      virtual ~BsAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob( const edm::EventSetup& ) ;
      virtual void endJob() ;

   private:
   
     //
     std::string fOutputFileName ;
     std::string theSrc ;
     TFile*      fOutputFile ;
     TH1D*       hGeneralId ;
     TH1D*       hPtbs ;
     TH1D*       hPbs ;
     TH1D*       hPhibs ;
     TH1D*       hEtabs ;
     TH1D*       hmumuMass ;
     TH1D*       hIdBsDaugs ;
     TH1D*       hIdBDaugs ;
    
     ofstream*   decayed; 
     ofstream*   undecayed; 
     int         nevent, nbs;
     
};

#endif
