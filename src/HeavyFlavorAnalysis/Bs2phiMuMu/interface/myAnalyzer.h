// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "HepMC/GenVertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetfwd.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "PhysicsTools/CandAlgos/interface/CloneProducer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
 
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "TFile.h"
#include "TTree.h"
#include <TLorentzVector.h>

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class myAnalyzer : public edm::EDAnalyzer {
 public:
  explicit myAnalyzer(const edm::ParameterSet&);
  ~myAnalyzer();
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillGeneratorBlock(const edm::Event&);
  virtual void fillRecTracks(const edm::Event&);
  virtual void fillJets(const edm::Event&);
  virtual void fillTrigger(const edm::Event&);
  virtual void fillMuonBlock(const edm::Event&);
  virtual void fillCands(const edm::Event&);
  virtual void endJob() ;
  float isolation(const edm::Event &iEvent, const reco::Candidate* aCand, double cone = 1.0);

  edm::ESHandle<TransientTrackBuilder> theB;
  edm::TriggerNames fTriggerNames;
  int HltEvtCnt;	       

  std::string fGenCandidatesLabel, fGenEventLabel, fTracksLabel,  fVertexLabel, fAssociatorLabel, fTrackingParticlesLabel;
  edm::InputTag fMuonLabel;
  edm::InputTag fJetsLabel;
  edm::InputTag fGenJetsLabel;
  edm::InputTag fHLTLabel;
  edm::InputTag fBCandLabel;
  edm::InputTag fBVtxLabel;
  bool storeTheBest;
  std::string whichBest;
  std::string outfile;

  // std::vector<const reco::Track*> RecTracks; 

  //static const int MAXEVENT = 1000000;
  static const int MAXGEN = 2000;
  static const int MAXBGEN = 5;
  static const int MAXTRK = 1000;
  static const int MAXMUON = 20;
  static const int MAXJET = 100;
  static const int MAXBRECO = 20000;

  int event;

  TFile* theFile;
  TTree* theTree; // event-wise tree

  Float_t  m_PVx;
  Float_t  m_PVy;
  Float_t  m_PVz;

  Int_t    m_Eve;
  Float_t  m_PtHat;

  Int_t    m_nGen;
  Int_t    m_pdgId[MAXGEN];
  Int_t    m_Number[MAXGEN];
  Int_t    m_Status[MAXGEN];
  Float_t  m_mcPx[MAXGEN];
  Float_t  m_mcPy[MAXGEN];
  Float_t  m_mcPz[MAXGEN];
  Float_t  m_mcE[MAXGEN];
  Float_t  m_mcPt[MAXGEN];
  Float_t  m_mcEta[MAXGEN];
  Float_t  m_mcPhi[MAXGEN];
  Float_t  m_mcVx[MAXGEN];
  Float_t  m_mcVy[MAXGEN];
  Float_t  m_mcVz[MAXGEN];

  Int_t    m_nbGen;
  Int_t    m_nBs;
  Float_t  m_bsGenMass[MAXBGEN];
  Float_t  m_phiGenMass[MAXBGEN];
  Float_t  m_mumuGenMass[MAXBGEN];

  Int_t    m_nTrk;
  Int_t    m_trkIndex[MAXTRK];
  Int_t    m_trkCharge[MAXTRK];
  Float_t  m_trkDof[MAXTRK];
  Int_t    m_trkNHit[MAXTRK];
  Float_t  m_trkPt[MAXTRK];
  Float_t  m_trkEta[MAXTRK];
  Float_t  m_trkPhi[MAXTRK];
  Float_t  m_trkIPt[MAXTRK];
  Float_t  m_trkIPtErr[MAXTRK];
  Float_t  m_trkIPz[MAXTRK];
  Float_t  m_trkIPzErr[MAXTRK];
  Float_t  m_trkChi2[MAXTRK];

  Int_t    m_nMu;
  Int_t    m_muIndex[MAXMUON];
  Float_t  m_muPt[MAXMUON];
  Float_t  m_muEta[MAXMUON];
  Float_t  m_muPhi[MAXMUON];
  Float_t  m_muIPt[MAXMUON];
  Float_t  m_muIPz[MAXMUON];
  Float_t  m_muChi2[MAXMUON];
  Int_t    m_muCharge[MAXMUON];
  Float_t  m_muDof[MAXMUON];
  Int_t    m_muNHit[MAXMUON];

  Int_t    m_nJet;
  Float_t  m_jtEne[MAXJET];
  Float_t  m_jtPhi[MAXJET];
  Float_t  m_jtEta[MAXJET];

  Int_t    m_nGJet;
  Float_t  m_gjtEne[MAXJET];
  Float_t  m_gjtPhi[MAXJET];
  Float_t  m_gjtEta[MAXJET];

  Int_t    *trigflag;
  
  Int_t    m_bsCand;
  Float_t  m_mumuMass[MAXBRECO];
  Float_t  m_phiMass[MAXBRECO];
  Float_t  m_bsMass[MAXBRECO];
  Float_t  m_bsIsoDR12[MAXBRECO];
  Float_t  m_bsIsoDR1[MAXBRECO];
  Float_t  m_bsIsoDR07[MAXBRECO];
  Float_t  m_bsIsoDR05[MAXBRECO];
  Int_t    m_bsVtx;
  Float_t  m_mumuDRVtx[MAXBRECO];
  Float_t  m_mumuMassVtx[MAXBRECO];
  Float_t  m_phiMassVtx[MAXBRECO];
  Float_t  m_bsMassVtx[MAXBRECO];
  Float_t  m_bsPtVtx[MAXBRECO];
  Float_t  m_bsPhiVtx[MAXBRECO];
  Float_t  m_bsEtaVtx[MAXBRECO];
  Float_t  m_chi2Vtx[MAXBRECO];
  Float_t  m_probVtx[MAXBRECO];
  Float_t  m_SVx[MAXBRECO];
  Float_t  m_SVy[MAXBRECO];
  Float_t  m_SVz[MAXBRECO];
  Float_t  m_decLen2D[MAXBRECO];
  Float_t  m_decSigma2D[MAXBRECO];
  Float_t  m_decLen3D[MAXBRECO];
  Float_t  m_decSigma3D[MAXBRECO];
  Float_t  m_bsPointing[MAXBRECO];
  Float_t  m_bsCosPointing[MAXBRECO];
};

