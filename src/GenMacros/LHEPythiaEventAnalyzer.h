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
     int         whichWeight;
     TFile*      fOutputFile ;
     // bool        isVBF;
     // int         lookFor;

     TH1D*       mZ;
     TH1D*       pTZ ;	   
     TH1D*       nJets ;
     TH1D*       pTjet1 ;
     TH1D*       pTjet2 ;
     TH1D*       pTjet3 ;
     TH1D*       pTZ_lowmass ;	   
     TH1D*       nJets_lowmass ;
     TH1D*       pTjet1_lowmass ;
     TH1D*       pTjet2_lowmass ;
     TH1D*       pTjet3_lowmass ;
    
     GreaterByPt<reco::GenJet> pTComparator_;
     int         nevent, neventpt0 ;
     edm::EDGetTokenT< LHEEventProduct > lhep_token;
     edm::EDGetTokenT<std::vector<reco::GenParticle> > genp_token;
   
};

#endif
