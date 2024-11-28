#ifndef LHEPythiaEventAnalyzer_H
#define LHEPythiaEventAnalyzer_H

// #include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// forward declarations
class TFile;
class TH1D;

class LHEPythiaEventAnalyzer : public edm::one::EDAnalyzer<> {

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
     int         whichWeight;
     TFile*      fOutputFile ;
     // bool        isVBF;
     // int         lookFor;

     TH1D*       helAngleKPlus;
    
     GreaterByPt<reco::GenJet> pTComparator_;
     int         nevent, neventpt0 ;
     edm::EDGetTokenT< LHEEventProduct > lhep_token;
     edm::EDGetTokenT<std::vector<reco::GenParticle> > genp_token;
   
};

#endif
