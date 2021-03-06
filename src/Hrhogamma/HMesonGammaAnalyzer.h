#ifndef HMesonGammaAnalyzer_H
#define HMesonGammaAnalyzer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// forward declarations
class TFile;
class TH1D;
class TEfficiency;

class HMesonGammaAnalyzer : public edm::EDAnalyzer
{

   public:
   
      //
      explicit HMesonGammaAnalyzer( const edm::ParameterSet& ) ;
      virtual ~HMesonGammaAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endJob() ;
      

   private:
   
      double deltaR(const double, const double, const double, const double);
      float ea (int, float);
      bool trackCuts(pat::PackedCandidate);
      std::string fOutputFileName ;
      edm::InputTag thePhotSrc ;
      edm::InputTag theTrkSrc ;
      edm::InputTag theTrkSrc2 ;
      edm::InputTag trigSrc ;
      std::string thePath ;

      //  
      TH1D*       MTktk ;
      TH1D*       Minv ;
      TH1D*       Minv_passHLT ;
      TH1D*       pTgamma ;
      TEfficiency* effpTgamma ;
      TH1D*       pTmeson ;
      TEfficiency* effpTmeson ;
      TH1D*       etagamma ;
      TEfficiency* effetagamma ;
      TH1D*       etameson ;
      TEfficiency* effetameson ;
      
      GreaterByPt<pat::Photon> pTComparator_;
      
      edm::EDGetTokenT< std::vector<pat::PackedCandidate> > trk_token;
      edm::EDGetTokenT< std::vector<pat::PackedCandidate> > trk_token2;
      edm::EDGetToken phot_token;
      edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<edm::TriggerResults> trigToken_;
      edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
      Float_t rho_; // the rho variable
      double hadronMass; // pi or K

      TFile* fOutputFile ;
      const edm::TriggerNames* triggerNames;

      int nevents, neventsGammaPass, neventsTauPass, neventsGammaMatch, neventsTauMatch;

};

#endif
