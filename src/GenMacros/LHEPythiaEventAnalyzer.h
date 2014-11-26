#ifndef LHEPythiaEventAnalyzer_H
#define LHEPythiaEventAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// forward declarations
class TFile;
class TH1D;

class LHEPythiaEventAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit LHEPythiaEventAnalyzer( const edm::ParameterSet& ) ;
      virtual ~LHEPythiaEventAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endJob() ;

   private:
   
     //
     std::string fOutputFileName ;
     std::string theSrc ;
     TFile*      fOutputFile ;
     // bool        isVBF;
     // int         lookFor;

     TH1D*       pTZ ;	   
     TH1D*       nJets ;
     TH1D*       pTjet1 ;
     TH1D*       pTjet2 ;
     TH1D*       pTjet3 ;
    
     GreaterByPt<reco::GenJet> pTComparator_;
     int         nevent, neventpt0 ;
   
};

#endif
