#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/Exception.h"

//#include "DQMOffline/RecoB/interface/JetTagPlotter.h"
//#include "DQMOffline/RecoB/interface/TagInfoPlotterFactory.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "DataFormats/Common/interface/View.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Validation/RecoB/plugins/myBTagAnalyzerMC.h"

using namespace reco;
using namespace edm;
using namespace std;
using namespace RecoBTag;
//using namespace BTagMCTools;

typedef std::pair<Jet, reco::JetFlavour> JetWithFlavour;

myBTagAnalyzerMC::myBTagAnalyzerMC(const edm::ParameterSet& pSet) :
  jetSelector(
    pSet.getParameter<double>("etaMin"),
    pSet.getParameter<double>("etaMax"),
    pSet.getParameter<double>("ptRecJetMin"),
    pSet.getParameter<double>("ptRecJetMax"),
    0.0, 99999.0,
    pSet.getParameter<double>("ratioMin"),
    pSet.getParameter<double>("ratioMax")
  ),
  jetMCSrc(pSet.getParameter<edm::InputTag>("jetMCSrc")),
  slInfoTag(pSet.getParameter<edm::InputTag>("softLeptonInfo")),
  moduleConfig(pSet.getParameter< edm::ParameterSet >("tagConfig")),
  jetCorrector(pSet.getParameter<std::string>("jetCorrection")),
  jetMatcher(pSet.getParameter<edm::ParameterSet>("recJetMatching"))
{
  vtxColl_ = pSet.getParameter<edm::InputTag>("vtxColl");
  jetColl_ = pSet.getParameter<edm::InputTag>("jetColl");
  electronColl_ = pSet.getParameter<edm::InputTag>("electronColl");
  muonColl_ = pSet.getParameter<edm::InputTag>("muonColl");
  metColl_ = pSet.getParameter<edm::InputTag>("metColl");
  // uncorrmetColl_ = pSet.getParameter<edm::InputTag>("uncorrmetColl");
  NeventsTOT_ = pSet.getParameter<int>( "NeventsTOT" );
  xsec_= pSet.getParameter<double>( "xsec" );
  lumi_= pSet.getParameter<double>( "lumi" );

  double ptRecJetMin = pSet.getParameter<double>("ptRecJetMin");
  jetMatcher.setThreshold(0.25 * ptRecJetMin);

  // True number of interaction for data produced as in: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
  TFile *da_=new TFile ("/afs/cern.ch/work/c/covarell/ttbar-ctag-ana/CMSSW_5_3_19/src/ExoDiBosonResonances/PATtupleProduction/src/MyDataPileupHistogram_True.root");
  TH1F *da = (TH1F*) da_->Get("pileup");
  
  // MC distribution of true number of interactions as in: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
  Double_t dat[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06 };
  
  // PileUp weights calculation
  double d,m;
  std::vector< float > mcNum;
  std::vector< float > dataNum;
  for (Int_t i=1; i< 50; i++){
    m=dat[i-1];
    d=da->GetBinContent(i);
    mcNum.push_back(m);
    dataNum.push_back(d);
  }
  LumiWeights_=edm::LumiReWeighting(mcNum, dataNum);
  
  eventInitialized = false;

  tagLabel = moduleConfig.getParameter<InputTag>("tagLabel");
  tagInfoLabel1 = moduleConfig.getParameter<InputTag>("tagInfoLabel1");
  tagInfoLabel2 = moduleConfig.getParameter<InputTag>("tagInfoLabel2");
  tagInfoLabel3 = moduleConfig.getParameter<InputTag>("tagInfoLabel3");

}

myBTagAnalyzerMC::~myBTagAnalyzerMC()
{  
}

void myBTagAnalyzerMC::beginJob()
{
  Service<TFileService> fs;

  theTree = fs->make<TTree>("theTree", "theTree");
  theTree->Branch("njet", &m_njet, "njet/i");
  theTree->Branch("njetinvmasses", &m_ninvmasses, "njetinvmasses/i");
  theTree->Branch("jetPt", m_jetPt, "jetPt[njet]/f");
  theTree->Branch("jetEta", m_jetEta, "jetEta[njet]/f");
  theTree->Branch("jetPhi", m_jetPhi, "jetPhi[njet]/f");
  theTree->Branch("jetMass", m_jetMass, "jetMass[njet]/f");
  theTree->Branch("jetDiscrCSV", m_jetDiscrCSV, "jetDiscrCSV[njet]/f");
  theTree->Branch("jetTrueFlavor", m_jetTrueFlavor, "jetTrueFlavor[njet]/f");
  theTree->Branch("jetNTrks", m_jetNTrks, "jetNTrks[njet]/i");
  theTree->Branch("jetNTrksSV", m_jetNTrksSV, "jetNTrksSV[njet]/i");
  theTree->Branch("jetWidth", m_jetWidth, "jetWidth[njet]/f");
  theTree->Branch("jetEccent", m_jetEccent, "jetEccent[njet]/f");
  theTree->Branch("jetMassSV", m_jetMassSV, "jetMassSV[njet]/f");
  theTree->Branch("jetAveTrkEtaRel", m_jetAveTrkEtaRel, "jetAveTrkEtaRel[njet]/f");
  theTree->Branch("jetESVOverE", m_jetESVOverE, "jetESVOverE[njet]/f");
  theTree->Branch("jetjetMass", m_jetjetMass, "jetjetMass[njetinvmasses]/f");
  theTree->Branch("jetjetMassIndex1", m_jetjetMassIndex1, "jetjetMassIndex1[njetinvmasses]/i");
  theTree->Branch("jetjetMassIndex2", m_jetjetMassIndex2, "jetjetMassIndex2[njetinvmasses]/i");
  theTree->Branch("dPhiLep1Met", &m_dPhiLep1Met, "dPhiLep1Met/f");
  theTree->Branch("dRLep1Met", &m_dRLep1Met, "dRLep1Met/f");
  theTree->Branch("lep1Pt", &m_lep1Pt, "lep1Pt/f");
  theTree->Branch("lep1Eta", &m_lep1Eta, "lep1Eta/f");
  theTree->Branch("lep1Phi", &m_lep1Phi, "lep1Phi/f");
  theTree->Branch("lep1Charge", &m_lep1Charge, "lep1Charge/f");
  theTree->Branch("lep1PFIso", &m_lep1PFIso, "lep1PFIso/f");
  theTree->Branch("lep1Type", &m_lep1Type, "lep1Type/f");
  theTree->Branch("met", &m_met, "met/f");
  theTree->Branch("metPhi", &m_metPhi, "metPhi/f");
  // theTree->Branch("uncorrmet", &m_uncorrmet, "uncorrmet/f");
  // theTree->Branch("uncorrmetPhi", &m_uncorrmetPhi, "uncorrmetPhi/f");
  theTree->Branch("NVertices", &m_NVertices, "NVertices/i");
  // theTree->Branch("PUWeight", &m_PUWeight, "PUWeight/f");
  theTree->Branch("metPx", &m_metPx, "metPx/f");
  theTree->Branch("metPy", &m_metPy, "metPy/f");
  theTree->Branch("NeventsTOT", &m_NeventsTOT, "NeventsTOT/i");
  theTree->Branch("xsec", &m_xsec, "xsec/d");
  theTree->Branch("lumi", &m_lumi, "lumi/d");
  theTree->Branch("weight", &m_weight, "weight/f");

}

unsigned int myBTagAnalyzerMC::SelectMuon(edm::Handle<reco::MuonCollection> muoH,
					  reco::Vertex primaryVertex, bool fully){
  unsigned int foundMu = 999;
  float ptMin = 0.;
  // for(reco::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
  for(unsigned imuo=0; imuo<muoH->size();++imuo) {

    reco::MuonRef muon(muoH,imuo);
    if(!(muon->isGlobalMuon())) continue;
    reco::TrackRef cktTrack = muon->track();
    if(cktTrack->pt()<20) continue;
    if(fabs(cktTrack->eta())>2.4) continue;
    if(fabs(cktTrack->phi())>3.2) continue;
    if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
    if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
    if(muon->numberOfMatches()<=1) continue;
    if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
    if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
    if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
    if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
    if(fully){if(MuonPFIso(muon, false)>0.2) continue;}
    if(muon->pt() > ptMin) {
      ptMin = muon->pt();
      foundMu = imuo;
    }
  }
  return foundMu;
}

unsigned int myBTagAnalyzerMC::SelectElectron(edm::Handle<reco::GsfElectronCollection> eleH,
			      reco::Vertex primaryVertex, 
			      edm::Handle<reco::ConversionCollection> hConversions, 
			      edm::Handle<reco::BeamSpot> beamSpotHandle,
			      const IsoDepositVals * electronIsoVals,
			      float rho, bool fully){
  unsigned int foundEle = 999;
  bool passEle=false;
  float ptMin = 0.;
 
  // for(reco::GsfElectronCollection::const_iterator electron = eleH->begin(); electron != eleH->end(); ++electron) {

  for(unsigned iele=0; iele<eleH->size();++iele) {
    reco::GsfElectronRef electron(eleH,iele);
    const reco::GsfElectron& anElectron = *electron;
    if(electron->pt()<20) continue;
    if(fully){
      if(ElectronPFIso(electron,electronIsoVals,rho)>0.1) continue;
      if(electron->pt()<20){if(ElectronPFIso(electron,electronIsoVals,rho)>0.07) continue;}
    }
    reco::GsfTrackRef gsfTrack = electron->gsfTrack();
    reco::SuperClusterRef superCluster = electron->superCluster();
    if(fabs(superCluster->eta())<=1.479){
      if(fabs(electron->deltaEtaSuperClusterTrackAtVtx())>=0.004) continue;
      if(fabs(electron->deltaPhiSuperClusterTrackAtVtx())>=0.030) continue;
      if(electron->sigmaIetaIeta()>=0.01) continue;
      if(electron->hadronicOverEm()>=0.12) continue;
      if(fabs(gsfTrack->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(gsfTrack->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(ConversionTools::hasMatchedConversion(anElectron, hConversions, beamSpotHandle->position())) continue;
      if(gsfTrack->trackerExpectedHitsInner().numberOfHits()!=0) continue;
      passEle=true;
    }
    if(fabs(superCluster->eta())>1.479 && fabs(superCluster->eta())<2.5){
      if(fabs(electron->deltaEtaSuperClusterTrackAtVtx())>=0.005) continue;
      if(fabs(electron->deltaPhiSuperClusterTrackAtVtx())>=0.020) continue;
      if(electron->sigmaIetaIeta()>=0.03) continue;
      if(electron->hadronicOverEm()>=0.10) continue;
      if(fabs(gsfTrack->dxy(primaryVertex.position()))>=0.02) continue;
      if(fabs(gsfTrack->dz(primaryVertex.position()))>=0.1) continue;
      if((fabs(1/electron->ecalEnergy() - electron->eSuperClusterOverP()/electron->ecalEnergy()))>=0.05) continue;
      if(ConversionTools::hasMatchedConversion(anElectron, hConversions, beamSpotHandle->position())) continue;
      if(gsfTrack->trackerExpectedHitsInner().numberOfHits()!=0) continue;
      passEle=true;
    }
    if(passEle && electron->pt() > ptMin) {
      foundEle = iele;
      ptMin = electron->pt();
    }
  }
  return foundEle;
}

float myBTagAnalyzerMC::MuonPFIso(reco::MuonRef muon, bool highpt){
  float sumChargedHadronPt = muon->pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon->pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon->pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon->pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/muon->pt();
  if(highpt){
    reco::TrackRef cktTrack = muon->track();
    iso = (sumChargedHadronPt+ max(0.,sumNeutralHadronEt+sumPhotonEt-0.5*sumPUPt))/cktTrack->pt();
  }
  return iso;
}


float myBTagAnalyzerMC::ElectronPFIso(reco::GsfElectronRef electron, const IsoDepositVals * electronIsoVals, float rho){

  float chargedHadronIso =  (*(*electronIsoVals)[0])[electron];
  float neutralHadronIso = (*(*electronIsoVals)[2])[electron];
  float photonIso =  (*(*electronIsoVals)[1])[electron];
  float thiseta = fabs(electron->superCluster()->eta());
  float Aeff=0.;
  if(thiseta<1.0) Aeff=0.13;
  if(thiseta>=1.0 && thiseta<1.479) Aeff=0.14;
  if(thiseta>=1.479 && thiseta<2.0) Aeff=0.07;
  if(thiseta>=2.0 && thiseta<2.2) Aeff=0.09;
  if(thiseta>=2.2 && thiseta<2.3) Aeff=0.11;
  if(thiseta>=2.3 && thiseta<2.4) Aeff=0.11;
  if(thiseta>=2.4) Aeff=0.14;
  float zero = 0.;
  float iso = (chargedHadronIso + max(zero, neutralHadronIso + photonIso - rho*Aeff))/electron->pt();
  return iso;
}

float myBTagAnalyzerMC::trackMom(float pt, float eta) {
  float theta = 2*atan(exp(-eta));
  return pt/sin(theta);
}

void myBTagAnalyzerMC::analyze(const edm::Event& iEvent, 
			       const edm::EventSetup& iSetup)
{

  int njet = 0;
  int ninvmasses = 0;

  for (int i = 0; i < NJETSMAX; i++) {
    m_jetPt[i] = -999.;
    m_jetEta[i] = -999.;
    m_jetPhi[i] = -999.;
    m_jetMass[i] = -999.;
    m_jetNTrks[i] = -999;
    m_jetNTrksSV[i] = -999;
    m_jetMassSV[i] = -999.;
    m_jetTrueFlavor[i] = -999.;
    m_jetDiscrCSV[i] = -999.;
    m_jetWidth[i] = -999.;
    m_jetEccent[i] = -999.;
    m_jetAveTrkEtaRel[i] = -999.;
    m_jetESVOverE[i] = -999.;
  }
  for (int j = 0; j < NJETINVMASSESMAX; j++) {
    m_jetjetMass[j] = -999.;
    m_jetjetMassIndex1[j] = -999;
    m_jetjetMassIndex2[j] = -999;
  }
  m_dPhiLep1Met= -999.;
  m_dRLep1Met= -999.;
  m_lep1Pt= -999.;
  m_lep1Phi= -999.;
  m_lep1Eta= -999.;
  m_lep1Charge= -999.;
  m_lep1PFIso= -999.;
  m_lep1Type= -999.;
  m_met= -999.;
  m_metPhi= -999.;
  // m_uncorrmet= -999.;
  // m_uncorrmetPhi= -999.;
  m_njet= 999;
  m_NVertices= 999;
  m_metPx= -999.;
  m_metPy= -999.;

  edm::Handle<vector<reco::GenParticle> > genParts;
  iEvent.getByLabel("genParticles", genParts);

  // *** Select lepton+jets only 
  bool isLeptonPlusJets = false;
  for(size_t ngenPart=0; ngenPart<genParts->size(); ngenPart++){
    const reco::GenParticle & genPart = (*genParts)[ngenPart];
    if ( (abs(genPart.pdgId()) == 13 || abs(genPart.pdgId()) == 11) && abs(genPart.mother()->pdgId()) == 24) {
      isLeptonPlusJets = true;
      break;
    }
  }
  if (!isLeptonPlusJets) return;
  // ***

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle); 

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vtxColl_, vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  edm::Handle<reco::PFJetCollection> ak5jetCands;
  iEvent.getByLabel(ak5JetColl_,ak5jetCands);

  edm::Handle<reco::MuonCollection> muoH;
  iEvent.getByLabel(muonColl_, muoH);

  edm::Handle<reco::GsfElectronCollection> eleH;
  iEvent.getByLabel(electronColl_, eleH);

  edm::Handle<std::vector<reco::PFMET> > met;
  iEvent.getByLabel(metColl_, met);

  //edm::Handle<std::vector<reco::PFMET> uncorrmet;
  // iEvent.getByLabel(uncorrmetColl_, uncorrmet);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel("kt6PFJets", "rho", rhoHandle);
  float rho = *(rhoHandle.product());

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);

  //Elementary method to get electron PFisolation ... mannaggia San Giorgino...
  std::vector<edm::InputTag> isovalvec;
  isovalvec.push_back(edm::InputTag("elPFIsoValueCharged03PFIdPFIso"));
  isovalvec.push_back(edm::InputTag("elPFIsoValueGamma03PFIdPFIso"));
  isovalvec.push_back(edm::InputTag("elPFIsoValueNeutral03PFIdPFIso"));    
  unsigned nTypes=3;
  IsoDepositVals electronIsoValPFId(nTypes);
  const IsoDepositVals * electronIsoVals = &electronIsoValPFId;
  for (size_t j = 0; j<isovalvec.size(); ++j) {
    iEvent.getByLabel(isovalvec[j], electronIsoValPFId[j]);
  }

  edm::Handle<JetFlavourMatchingCollection> jetMC;
  FlavourMap flavours;
  LeptonMap leptons;

  iEvent.getByLabel(jetMCSrc, jetMC);
  for (JetFlavourMatchingCollection::const_iterator iter = jetMC->begin();
       iter != jetMC->end(); ++iter) {
    unsigned int fl = std::abs(iter->second.getFlavour());
    flavours.insert(std::make_pair(iter->first, fl));
    const reco::JetFlavour::Leptons &lep = iter->second.getLeptons();
    leptons.insert(std::make_pair(iter->first, lep));
  }

  edm::Handle<reco::SoftLeptonTagInfoCollection> infoHandle;
  iEvent.getByLabel(slInfoTag, infoHandle);

  math::PtEtaPhiELorentzVector lep1;
  float lep1Type=0;
  float lep1Charge=0;
  float lep1PFIso=100.; 

  unsigned SelectedEle = SelectElectron(eleH, primaryVertex, hConversions, beamSpotHandle, electronIsoVals, rho, true);

  if (SelectedEle < 999) {
    reco::GsfElectronRef SelEle(eleH,SelectedEle);
    lep1 = SelEle->p4();
    lep1Charge=SelEle->charge();
    lep1PFIso=ElectronPFIso(SelEle,electronIsoVals,rho);
    lep1Type=1;
  }

  unsigned SelectedMuo = SelectMuon(muoH, primaryVertex, true);

  if (SelectedMuo < 999) {
    reco::MuonRef SelMuo(muoH,SelectedMuo); 
    lep1 = SelMuo->p4();
    lep1Charge=SelMuo->charge();
    lep1PFIso=MuonPFIso(SelMuo,true);
    lep1Type=2;
  }

  // cout << "Oppure quaggiu'? " << endl; 

  // Look first at the jetTags

  edm::Handle<reco::JetTagCollection> tagHandle;
  iEvent.getByLabel(tagLabel, tagHandle);
  const reco::JetTagCollection & tagColl = *(tagHandle.product());
  LogDebug("Info") << "Found " << tagColl.size() << " B candidates in collection " << tagLabel;

  for (JetTagCollection::const_iterator tagI = tagColl.begin();
       tagI != tagColl.end(); ++tagI) {
    // Identify parton associated to jet.
    
    JetWithFlavour jetWithFlavour;
    if (!getJetWithFlavour(tagI->first, flavours, jetWithFlavour, iSetup))
      continue;
    if (!jetSelector(jetWithFlavour.first, std::abs(jetWithFlavour.second.getFlavour()), infoHandle))
      continue;
    
    m_jetPt[njet]=jetWithFlavour.first.pt();
    m_jetEta[njet]=jetWithFlavour.first.eta();
    m_jetPhi[njet]=jetWithFlavour.first.phi();
    // cout << "loop on JetTag *** " << jetWithFlavour.first.pt() << " " << jetWithFlavour.first.eta() << " " << jetWithFlavour.first.phi() << endl; 
    m_jetMass[njet]=jetWithFlavour.first.mass();
    m_jetDiscrCSV[njet] = tagI->second;
    m_jetTrueFlavor[njet] = jetWithFlavour.second.getFlavour();
    
    int njet2 = 0;
    for (JetTagCollection::const_iterator tagII = tagColl.begin();
	 tagII != tagI; ++tagII) {
      
      JetWithFlavour jetWithFlavourII;
      if (!getJetWithFlavour(tagII->first, flavours, jetWithFlavourII, iSetup))
	continue;
      if (!jetSelector(jetWithFlavourII.first, std::abs(jetWithFlavourII.second.getFlavour()), infoHandle))
	continue;
      
      TLorentzVector lvjet1;
      lvjet1.SetPxPyPzE(jetWithFlavour.first.px(),jetWithFlavour.first.py(),
			jetWithFlavour.first.pz(),jetWithFlavour.first.energy() );
      TLorentzVector lvjet2;
      lvjet2.SetPxPyPzE(jetWithFlavourII.first.px(),jetWithFlavourII.first.py(),
			jetWithFlavourII.first.pz(),jetWithFlavourII.first.energy() );
      TLorentzVector jetjet = lvjet1 + lvjet2;
      if (jetjet.M() > 40 && jetjet.M() < 150) {
	m_jetjetMass[ninvmasses] = jetjet.M();
	m_jetjetMassIndex1[ninvmasses] = njet2;
	m_jetjetMassIndex2[ninvmasses] = njet;
	ninvmasses++;
      }
      njet2++;
    }
    njet++;
  }	

  // Now look at the TagInfos

  unsigned int nInputTags = 2;
  vector< edm::Handle< View<BaseTagInfo> > > tagInfoHandles;
  edm::Handle< View<BaseTagInfo> > tagInfoHandle;
  iEvent.getByLabel(tagInfoLabel1, tagInfoHandle);
  tagInfoHandles.push_back(tagInfoHandle);
  edm::Handle< View<BaseTagInfo> > tagInfoHandle2;
  iEvent.getByLabel(tagInfoLabel2, tagInfoHandle2);
  tagInfoHandles.push_back(tagInfoHandle2);
  edm::Handle< View<BaseTagInfo> > tagInfoHandle3;
  iEvent.getByLabel(tagInfoLabel3, tagInfoHandle3);
  // tagInfoHandles.push_back(tagInfoHandle3);
  
  edm::RefToBase<Jet> jetRef;
  for (unsigned int iTagInfo = 0; iTagInfo < nInputTags; iTagInfo++) {
    for (unsigned int ijetTagInfo = 0; ijetTagInfo < tagInfoHandles[iTagInfo]->size(); ijetTagInfo++) {
      // cout << "*** loop on TagInfo index " << iTagInfo << " " << ijetTagInfo << " " << tagInfoHandles[iTagInfo]->size() << endl;  
      const BaseTagInfo &baseTagInfo = (*tagInfoHandles[iTagInfo])[ijetTagInfo];
      jetRef = baseTagInfo.jet();
    
    // Identify parton associated to jet.
    
      JetWithFlavour jetWithFlavour;
      if (!getJetWithFlavour(jetRef, flavours, jetWithFlavour, iSetup))
	continue;
      if (!jetSelector(jetWithFlavour.first, std::abs(jetWithFlavour.second.getFlavour()), infoHandle))
	continue;
      
      // cout << "*** loop on TagInfo " << jetWithFlavour.first.pt() << " " << jetWithFlavour.first.eta() << " " << jetWithFlavour.first.phi() << endl;   
      
      // now find out which jet it is in the previous collection
      int theJet = -1;
      for (int ijet = 0; ijet < njet; ijet++) {
	if(fabs(m_jetEta[ijet]-jetWithFlavour.first.eta()) < 0.01 &&
	   fabs(m_jetPhi[ijet]-jetWithFlavour.first.phi()) < 0.01 &&
	   fabs(m_jetPt[ijet]-jetWithFlavour.first.pt()) < 0.1 ) theJet = ijet;
      }
      
      // cout << "*** loop on TagInfo " << theJet << endl;

      // const JetTagComputer::TagInfoHelper helper(baseTagInfo);
      // const TaggingVariableList& vars = computer->taggingVariables(helper);

      const reco::TaggingVariableList &vars = baseTagInfo.taggingVariables();
      if (vars.checkTag(getTaggingVariableName("vertexNTracks"))) 
	m_jetNTrksSV[theJet] = vars.get(getTaggingVariableName("vertexNTracks"));
      if (vars.checkTag(getTaggingVariableName("trackEtaRel"))) {
	std::vector<float> theTrackRaps = vars.getList(getTaggingVariableName("trackEtaRel"), false);
	float ave = 0.;
	for (unsigned int itr = 0; itr < theTrackRaps.size(); itr++) {
	  ave += theTrackRaps[itr];
	}
	ave = ave / theTrackRaps.size();
	m_jetAveTrkEtaRel[theJet] = ave;
      }

      float jetw = 0.;
      float jetecc = 0.;
      float sumpt = 0.;
      const reco::TrackIPTagInfo * tptagInfo = dynamic_cast<const reco::TrackIPTagInfo *>(&baseTagInfo);
      if (tptagInfo) {
	 std::vector<std::size_t> sortedIndices = tptagInfo->sortedIndexes(reco::TrackIPTagInfo::IP2DSig);
	 reco::TrackRefVector sortedTracks = tptagInfo->sortedTracks(sortedIndices);
	 m_jetNTrks[theJet] = sortedIndices.size();
	 for(unsigned int n=0; n != sortedIndices.size(); ++n) {
	    const reco::TrackRef& track = sortedTracks[n];
	    math::PtEtaPhiELorentzVector trackp4 = math::PtEtaPhiELorentzVector(track->pt(),track->eta(),track->phi(),trackMom(track->pt(),track->eta()) );
	    jetw += ROOT::Math::VectorUtil::DeltaR(jetWithFlavour.first.p4(),trackp4)*track->pt();
	    jetecc += pow(ROOT::Math::VectorUtil::DeltaR(jetWithFlavour.first.p4(),trackp4),2)*track->pt();
	    sumpt += track->pt();
	 }
         m_jetWidth[theJet] = jetw/sumpt;
         m_jetEccent[theJet] = sqrt(jetecc/sumpt);  
      }

      math::PtEtaPhiELorentzVector alltracks = math::PtEtaPhiELorentzVector(0.,0.,0.,0.);
      float energyalltracks = 0.;
      const reco::SecondaryVertexTagInfo * svtagInfo = dynamic_cast<const reco::SecondaryVertexTagInfo *>(&baseTagInfo);
      if (svtagInfo) {
	 const reco::TrackRefVector vertexTracks = svtagInfo->vertexTracks();
	 for(unsigned int n=0; n != vertexTracks.size(); ++n) {
	    const reco::TrackRef& track = vertexTracks[n];
	    math::PtEtaPhiELorentzVector trackp4 = math::PtEtaPhiELorentzVector(track->pt(),track->eta(),track->phi(),trackMom(track->pt(),track->eta()) );
	    alltracks += trackp4;
	    energyalltracks += trackMom(track->pt(),track->eta()) ;
	 }
         m_jetMassSV[theJet] = alltracks.M();
         m_jetESVOverE[theJet] = energyalltracks/jetWithFlavour.first.energy();  
      } 	    
    } 
  }

  TLorentzVector MET;
  MET.SetPxPyPzE(met->begin()->px(),met->begin()->py(),met->begin()->pz(),met->begin()->energy());

  m_lep1Pt=lep1.pt();
  m_lep1Eta=lep1.eta();
  m_lep1Phi=lep1.phi();
  m_lep1Charge=lep1Charge;
  m_lep1PFIso=lep1PFIso;
  m_lep1Type=lep1Type;
  m_dPhiLep1Met=ROOT::Math::VectorUtil::DeltaPhi(met->begin()->p4(),lep1);
  m_dRLep1Met=ROOT::Math::VectorUtil::DeltaR(met->begin()->p4(),lep1);
  m_met=MET.Pt();
  m_metPhi=met->begin()->phi();
  // m_uncorrmet=uncorrmet->begin()->pt();
  // m_uncorrmetPhi=uncorrmet->begin()->phi();
  m_njet=njet;
  m_ninvmasses=ninvmasses;
  m_NVertices=vertices->size();
  // m_PUWeight=MyWeight;
  m_metPx=MET.Px();
  m_metPy=MET.Py();
  m_NeventsTOT=NeventsTOT_;
  m_xsec=xsec_;
  m_lumi=lumi_;
  m_weight=xsec_*lumi_/NeventsTOT_;
  theTree->Fill();
  
}

bool myBTagAnalyzerMC::getJetWithFlavour( edm::RefToBase<Jet> jetRef, 
					   FlavourMap flavours,
					   JetWithFlavour & jetWithFlavour, 
					   const edm::EventSetup & es)
{

  edm::ProductID recProdId = jetRef.id();
  edm::ProductID refProdId = (flavours.begin() == flavours.end())
    ? recProdId
    : flavours.begin()->first.id();

  if (!eventInitialized) {
    jetCorrector.setEventSetup(es);
    if (recProdId != refProdId) {
      edm::RefToBaseVector<Jet> refJets;
      for(FlavourMap::const_iterator iter = flavours.begin();
          iter != flavours.end(); ++iter)
        refJets.push_back(iter->first);
      const edm::RefToBaseProd<Jet> recJetsProd(jetRef);
      edm::RefToBaseVector<Jet> recJets;
      for(unsigned int i = 0; i < recJetsProd->size(); i++)
        recJets.push_back(edm::RefToBase<Jet>(recJetsProd, i));
      jetMatcher.matchCollections(refJets, recJets, es);
    }
    eventInitialized = true;
  }

  if (recProdId != refProdId) {
    jetRef = jetMatcher(jetRef);
    if (jetRef.isNull())
      return false;
  }

  jetWithFlavour.first = jetCorrector(*jetRef);

  jetWithFlavour.second = reco::JetFlavour(jetWithFlavour.first.p4(), math::XYZPoint (0,0,0), flavours[jetRef]);

  LogTrace("Info") << "Found jet with flavour "<<jetWithFlavour.second.getFlavour()<<endl;
  LogTrace("Info") << jetWithFlavour.first.p()<<" , "<< jetWithFlavour.first.pt()<<" - "
   << jetWithFlavour.second.getLorentzVector().P()<<" , "<< jetWithFlavour.second.getLorentzVector().Pt()<<endl;

  return true;
}

void myBTagAnalyzerMC::endJob()
{
}


//define this as a plug-in
DEFINE_FWK_MODULE(myBTagAnalyzerMC);
