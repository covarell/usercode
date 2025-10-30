#ifndef HPhiPhiAnalyzer_H
#define HPhiPhiAnalyzer_H

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// forward declarations
class TFile;
class TH1D;
class TEfficiency;

class HPhiPhiAnalyzer : public edm::one::EDAnalyzer<> {

   public:
   
      //
      explicit HPhiPhiAnalyzer( const edm::ParameterSet& ) ;
      virtual ~HPhiPhiAnalyzer() {} // no need to delete ROOT stuff
                               // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endJob() ;
      

   private:
   
      double deltaR(const double, const double, const double, const double);
      bool trackCuts(pat::PackedCandidate);
      std::string fOutputFileName ;
      edm::InputTag theTrkSrc ;
      edm::InputTag theTrkSrc2 ;
      edm::InputTag trigSrc ;
      std::string thePath ;

      //  
      TH1D*       MTktk1 ;
      TH1D*       MTktk2 ;
      TH1D*       Minv ;
      TH1D*       Minv_passHLT ;
      TH1D*       pTmeson ;
      TEfficiency* effpTmeson ;
      TH1D*       etameson ;
      TEfficiency* effetameson ;
      
            
      edm::EDGetTokenT< std::vector<pat::PackedCandidate> > trk_token;
      edm::EDGetTokenT< std::vector<pat::PackedCandidate> > trk_token2;
      edm::EDGetTokenT<edm::TriggerResults> trigToken_;
      edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
      
      double hadronMass; // pi or K

      TFile* fOutputFile ;
      const edm::TriggerNames* triggerNames;

  int nevents, neventsTauPass, neventsTauMatch;

};

#endif
