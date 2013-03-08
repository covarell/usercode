#ifndef LHEEventAnalyzer_H
#define LHEEventAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

// forward declarations
class TFile;
class TH2D;

class LHEEventAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit LHEEventAnalyzer( const edm::ParameterSet& ) ;
      virtual ~LHEEventAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endJob() ;

   private:
   
     //
     std::string fOutputFileName ;
     std::string theSrc ;
     TFile*      fOutputFile ;
     TH2D*       Pt_sigVBF ;	   
     TH2D*       Y_sigVBF ;
     TH2D*       Pt_bkg ;
     TH2D*       PtOverM_bkg ;
     TH2D*       Y_bkg ;
    
     int         nevent, neventpt0 ;
     bool        isVBF;
     
};

#endif
