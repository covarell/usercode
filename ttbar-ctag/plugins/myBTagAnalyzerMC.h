#ifndef myBTagAnalyzerMC_H
#define myBTagAnalyzerMC_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// #include "DQMOffline/RecoB/interface/BaseBTagPlotter.h"
#include "DataFormats/Common/interface/ValueMap.h" 
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" 
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h" 
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
//#include "DQMOffline/RecoB/interface/BTagDifferentialPlot.h"
#include "DQMOffline/RecoB/interface/AcceptJet.h"
//#include "DQMOffline/RecoB/interface/JetTagPlotter.h"
//#include "DQMOffline/RecoB/interface/TagCorrelationPlotter.h"
//#include "DQMOffline/RecoB/interface/BaseTagInfoPlotter.h"
#include "DQMOffline/RecoB/interface/Tools.h"
//#include "RecoBTag/MCTools/interface/JetFlavourIdentifier.h"
//#include "RecoBTag/MCTools/interface/JetFlavour.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "DQMOffline/RecoB/interface/CorrectJet.h"
#include "DQMOffline/RecoB/interface/MatchJet.h"

#include "TTree.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>
#include <map>

//class CaloJetRef;

/** \class myBTagAnalyzerMC
 *
 *  Top level steering routine for b tag performance analysis.
 *
 */

class myBTagAnalyzerMC : public edm::EDAnalyzer {
   public:
      explicit myBTagAnalyzerMC(const edm::ParameterSet& pSet);

      ~myBTagAnalyzerMC();

      virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      virtual void beginJob() ;
      virtual void endJob();

   private:

  struct JetRefCompare :
       public std::binary_function<edm::RefToBase<reco::Jet>, edm::RefToBase<reco::Jet>, bool> {
    inline bool operator () (const edm::RefToBase<reco::Jet> &j1,
                             const edm::RefToBase<reco::Jet> &j2) const
    { return j1.id() < j2.id() || (j1.id() == j2.id() && j1.key() < j2.key()); }
  };

  // Get histogram plotting options from configuration.
  typedef std::pair<reco::Jet, reco::JetFlavour> JetWithFlavour;
  typedef std::map<edm::RefToBase<reco::Jet>, unsigned int, JetRefCompare> FlavourMap;
  typedef std::map<edm::RefToBase<reco::Jet>, reco::JetFlavour::Leptons, JetRefCompare> LeptonMap;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

  //  reco::JetFlavour getJetFlavour(
  //	edm::RefToBase<reco::Jet> caloRef, FlavourMap flavours);
  bool getJetWithFlavour( edm::RefToBase<reco::Jet> caloRef,
                         FlavourMap flavours, JetWithFlavour &jetWithFlavour,
                         const edm::EventSetup & es);
  unsigned int SelectMuon(edm::Handle<reco::MuonCollection> muoH,
		  reco::Vertex primaryVertex, bool fully);
  unsigned int SelectElectron(edm::Handle<reco::GsfElectronCollection> eleH,
		      reco::Vertex primaryVertex, edm::Handle<reco::ConversionCollection> hConversions, edm::Handle<reco::BeamSpot> beamSpotHandle, const IsoDepositVals * electronIsoVals, float rho, bool fully);
  float MuonPFIso(reco::MuonRef muon, bool highpt);
  float ElectronPFIso(reco::GsfElectronRef electron, const IsoDepositVals * electronIsoVals, float rho);
  float trackMom(float pt, float eta);

  AcceptJet jetSelector;   // Decides if jet and parton satisfy kinematic cuts.
  std::vector<double> etaRanges, ptRanges;
  edm::InputTag jetMCSrc;
  edm::InputTag slInfoTag;
  edm::ParameterSet moduleConfig;

  edm::InputTag tagLabel;
  edm::InputTag tagInfoLabel1;
  edm::InputTag tagInfoLabel2;
  edm::InputTag tagInfoLabel3;

  bool eventInitialized;
  CorrectJet jetCorrector;
  MatchJet jetMatcher;

  // Variables in tree
  static const int NJETSMAX = 20; 
  static const int NJETINVMASSESMAX = 40;

  float m_genEvent;

  TTree *theTree;
  float m_dPhiLep1Met;
  float m_dRLep1Met;
  float m_lep1Pt;
  float m_lep1Eta;
  float m_lep1Phi;
  float m_lep1Charge;
  float m_lep1PFIso;
  float m_lep1Type;
  float m_met;
  float m_metPhi;
  float m_uncorrmet;
  float m_uncorrmetPhi;
  int m_njet;
  int m_ninvmasses;
  float m_jetPt[NJETSMAX];
  float m_jetEta[NJETSMAX];
  float m_jetPhi[NJETSMAX];
  float m_jetMass[NJETSMAX];
  int m_jetNTrks[NJETSMAX];
  int m_jetNTrksSV[NJETSMAX];
  float m_jetMassSV[NJETSMAX];
  float m_jetTrueFlavor[NJETSMAX];
  float m_jetDiscrCSV[NJETSMAX];
  float m_jetWidth[NJETSMAX];
  float m_jetEccent[NJETSMAX];
  float m_jetAveTrkEtaRel[NJETSMAX];
  float m_jetESVOverE[NJETSMAX];
  float m_jetjetMass[NJETINVMASSESMAX];
  int m_jetjetMassIndex1[NJETINVMASSESMAX];
  int m_jetjetMassIndex2[NJETINVMASSESMAX];
  int m_NVertices;
  // float m_PUWeight;
  float m_metPx;
  float m_metPy;
  int m_NeventsTOT;
  double m_xsec;
  double m_lumi;
  double m_weight;

  edm::LumiReWeighting LumiWeights_;
  // bool isData;
  edm::InputTag vtxColl_;
  edm::InputTag jetColl_;
  edm::InputTag electronColl_;
  edm::InputTag muonColl_;
  edm::InputTag ak5JetColl_;  // useful?
  edm::InputTag metColl_;
  edm::InputTag uncorrmetColl_;
  int NeventsTOT_;
  double xsec_;
  double lumi_;
};


#endif
