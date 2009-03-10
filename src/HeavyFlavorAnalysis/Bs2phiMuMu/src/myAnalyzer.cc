
#include <iostream>
#include <TRandom.h>
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "HeavyFlavorAnalysis/Bs2phiMuMu/interface/myAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace HepMC;

float pai = 3.1415926;
int event = 0;
// ----------------------------------------------------------------------
myAnalyzer::myAnalyzer(const ParameterSet& iConfig):
  
  fGenEventLabel(iConfig.getUntrackedParameter<string>("generatorEvent", string("source"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("gsWithMaterialTracks"))),
  fVertexLabel(iConfig.getUntrackedParameter<string>("vertexLabel", string("theVertex"))), 
  fAssociatorLabel(iConfig.getUntrackedParameter<string>("associatorLabel", string("TrackAssociatorByChi2"))), 
  fTrackingParticlesLabel(iConfig.getUntrackedParameter<string>("trackingParticlesLabel", string("trackingtruthprod"))),
  fMuonLabel(iConfig.getParameter<InputTag>("muonLabel")),
  fJetsLabel(iConfig.getParameter<InputTag>("jetsLabel")),
  fGenJetsLabel(iConfig.getParameter<InputTag>("genjetsLabel")),
  fHLTLabel(iConfig.getParameter<InputTag>("HLTLabel")),
  fBCandLabel(iConfig.getParameter<InputTag>("bCandLabel")),
  fBVtxLabel(iConfig.getParameter<InputTag>("bVtxLabel")),
  storeTheBest(iConfig.getParameter<bool>("storeTheBest")),
  whichBest(iConfig.getUntrackedParameter<string>("whichBest",string("none"))),
  outfile(iConfig.getParameter<string>("outfile"))

{

  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "-------------- myAnalyzer constructor: " << endl;
  cout << "OUTPUT FILE        : " << outfile << endl;
  cout << "- GenEvent label   : " << fGenEventLabel << endl;
  cout << "- Tracks label     : " << fTracksLabel << endl;
  cout << "- Muons label      : " << fMuonLabel << endl;
  cout << "- GenJets label    : " << fGenJetsLabel << endl;
  cout << "- Jets label       : " << fJetsLabel << endl;
  cout << "- Vertex label     : " << fVertexLabel << endl;
  cout << "- Associator label : " << fAssociatorLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  // create ROOT files
  theFile = new TFile(outfile.c_str(),"recreate");
  theFile->cd();
  theTree  = new TTree("ntp1","bkg tree"); 

  // define branches and leaves:
  // primary vertex
  theTree->Branch("PVx", &m_PVx, "PVx/F");
  theTree->Branch("PVy", &m_PVy, "PVy/F");
  theTree->Branch("PVz", &m_PVz, "PVz/F");

  theTree->Branch("Eve", &m_Eve, "Eve/I");
  theTree->Branch("PtHat", &m_PtHat, "PtHat/F");

  theTree->Branch("nGen", &m_nGen, "nGen/I");
  theTree->Branch("pdgId", m_pdgId, "pdgId[nGen]/I");
  theTree->Branch("Number", m_Number, "Number[nGen]/I");
  theTree->Branch("Status", m_Status, "Status[nGen]/I");
  theTree->Branch("mcPx", m_mcPx, "mcPx[nGen]/F");
  theTree->Branch("mcPy", m_mcPy, "mcPy[nGen]/F");
  theTree->Branch("mcPz", m_mcPz, "mcPz[nGen]/F");
  theTree->Branch("mcE", m_mcE, "mcE[nGen]/F");
  theTree->Branch("mcPt", m_mcPt, "mcPt[nGen]/F");
  theTree->Branch("mcEta", m_mcEta, "mcEta[nGen]/F");
  theTree->Branch("mcPhi", m_mcPhi, "mcPhi[nGen]/F");
  theTree->Branch("mcVx", m_mcVx, "mcVx[nGen]/F");
  theTree->Branch("mcVy", m_mcVy, "mcVy[nGen]/F");
  theTree->Branch("mcVz", m_mcVz, "mcVz[nGen]/F");
  // generated bs
  theTree->Branch("nBs", &m_nBs, "nBs/I");
  theTree->Branch("nbGen", &m_nbGen, "nbGen/I");
  theTree->Branch("bsGenMass", m_bsGenMass, "bsGenMass[nbGen]/F");
  theTree->Branch("phiGenMass", m_phiGenMass, "phiGenMass[nbGen]/F");
  theTree->Branch("mumuGenMass", m_mumuGenMass, "mumuGenMass[nbGen]/F");
  //tracks
  theTree->Branch("nTrk", &m_nTrk, "nTrk/I");
  theTree->Branch("trkIndex", m_trkIndex, "trkIndex[nTrk]/I");
  theTree->Branch("trkPt", m_trkPt, "trkPt[nTrk]/F");
  theTree->Branch("trkEta", m_trkEta, "trkEta[nTrk]/F");
  theTree->Branch("trkPhi", m_trkPhi, "trkPhi[nTrk]/F");
  theTree->Branch("trkIPt", m_trkIPt, "trkIPt[nTrk]/F");
  theTree->Branch("trkIPtErr", m_trkIPtErr, "trkIPtErr[nTrk]/F");
  theTree->Branch("trkIPz", m_trkIPz, "trkIPz[nTrk]/F");
  theTree->Branch("trkIPzErr", m_trkIPzErr, "trkIPzErr[nTrk]/F");
  theTree->Branch("trkChi2", m_trkChi2, "trkChi2[nTrk]/F");
  theTree->Branch("trkCharge", m_trkCharge, "trkCharge[nTrk]/I");
  theTree->Branch("trkDof", m_trkDof, "trkDof[nTrk]/I");
  theTree->Branch("trkNHit", m_trkNHit, "trkNHit[nTrk]/I");
  // muons
  theTree->Branch("nMu", &m_nMu, "nMu/I");
  theTree->Branch("muIndex", m_muIndex, "muIndex[nMu]/I");
  theTree->Branch("muPt", m_muPt, "muPt[nMu]/F");
  theTree->Branch("muEta", m_muEta, "muEta[nMu]/F");
  theTree->Branch("muPhi", m_muPhi, "muPhi[nMu]/F");
  theTree->Branch("muIPt", m_muIPt, "muIPt[nMu]/F");
  theTree->Branch("muIPz", m_muIPz, "muIPz[nMu]/F");
  theTree->Branch("muChi2", m_muChi2, "muChi2[nMu]/F");
  theTree->Branch("muCharge", m_muCharge, "muCharge[nMu]/I");
  theTree->Branch("muDof", m_muDof, "muDof[nMu]/I");
  theTree->Branch("muNHit", m_muNHit, "muNHit[nMu]/I");
  // gen jets
  theTree->Branch("nGJet", &m_nGJet, "nGJet/I");
  theTree->Branch("gjtEne", m_gjtEne, "gjtEne[nGJet]/F");
  theTree->Branch("gjtEta", m_gjtEta, "gjtEta[nGJet]/F");
  theTree->Branch("gjtPhi", m_gjtPhi, "gjtPhi[nGJet]/F");
  // jets
  theTree->Branch("nJet", &m_nJet, "nJet/I");
  theTree->Branch("jtEne", m_jtEne, "jtEne[nJet]/F");
  theTree->Branch("jtEta", m_jtEta, "jtEta[nJet]/F");
  theTree->Branch("jtPhi", m_jtPhi, "jtPhi[nJet]/F");
  // HLT 
  HltEvtCnt = 0;
  static const int MAXTRG = 200;
  trigflag = new int[MAXTRG];

  //Bs bare candidates
  theTree->Branch("bsCand", &m_bsCand, "bsCand/I");
  theTree->Branch("mumuMass", m_mumuMass, "mumuMass[bsCand]/F");
  theTree->Branch("phiMass", m_phiMass, "phiMass[bsCand]/F");
  theTree->Branch("bsMass", m_bsMass, "bsMass[bsCand]/F");
     
  //Bs vertexed candidates  
  theTree->Branch("bsVtx", &m_bsVtx, "bsVtx/I");
  theTree->Branch("mumuMassVtx", m_mumuMassVtx, "mumuMassVtx[bsVtx]/F");
  theTree->Branch("mumuDRVtx", m_mumuDRVtx, "mumuDRVtx[bsVtx]/F");
  theTree->Branch("phiMassVtx", m_phiMassVtx, "phiMassVtx[bsVtx]/F");
  theTree->Branch("bsMassVtx", m_bsMassVtx, "bsMassVtx[bsVtx]/F");
  theTree->Branch("bsPtVtx", m_bsPtVtx, "bsPtVtx[bsVtx]/F");
  theTree->Branch("bsPhiVtx", m_bsPhiVtx, "bsPhiVtx[bsVtx]/F");
  theTree->Branch("bsEtaVtx", m_bsEtaVtx, "bsEtaVtx[bsVtx]/F");
  theTree->Branch("chi2Vtx", m_chi2Vtx, "chi2Vtx[bsVtx]/F");
  theTree->Branch("probVtx", m_probVtx, "probVtx[bsVtx]/F");
  theTree->Branch("SVx", m_SVx, "SVx[bsVtx]/F");
  theTree->Branch("SVy", m_SVy, "SVy[bsVtx]/F");
  theTree->Branch("SVz", m_SVz, "SVz[bsVtx]/F");
  theTree->Branch("decLen2D", m_decLen2D, "decLen2D[bsVtx]/F");
  theTree->Branch("decSigma2D", m_decSigma2D, "decSigma2D[bsVtx]/F");
  theTree->Branch("decLen3D", m_decLen3D, "decLen3D[bsVtx]/F");
  theTree->Branch("decSigma3D", m_decSigma3D, "decSigma3D[bsVtx]/F");
  theTree->Branch("bsPointing", m_bsPointing, "bsPointing[bsVtx]/F");
  theTree->Branch("bsCosPointing", m_bsCosPointing, "bsCosPointing[bsVtx]/F");
  theTree->Branch("bsIsoDR12", m_bsIsoDR12, "bsIsoDR12[bsVtx]/F");
  theTree->Branch("bsIsoDR1", m_bsIsoDR1, "bsIsoDR1[bsVtx]/F");
  theTree->Branch("bsIsoDR07", m_bsIsoDR07, "bsIsoDR07[bsVtx]/F");
  theTree->Branch("bsIsoDR05", m_bsIsoDR05, "bsIsoDR05[bsVtx]/F");
}

// ----------------------------------------------------------------------
myAnalyzer::~myAnalyzer() {  
}

// ----------------------------------------------------------------------
void myAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup) {

  // -- Get generator block directly
  fillGeneratorBlock(iEvent);
  
  // -- Get Tracks Block
  fillRecTracks(iEvent);

  // -- Get the Muon block
  fillMuonBlock(iEvent);

  // -- Get the Muon block
  fillJets(iEvent);
 
  // -- Get the Muon block
  fillTrigger(iEvent);

  // -- Form candidates 
  fillCands(iEvent);

  theTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void myAnalyzer::fillGeneratorBlock(const edm::Event &iEvent) {

  Handle<HepMCProduct> evt;

  double genEventScale;
  event++;
  m_Eve = event;
  
  try {

    iEvent.getByLabel(fGenEventLabel.c_str(), evt); 
    const HepMC::GenEvent *genEvent = evt->GetEvent();
  
    genEventScale = (*genEvent).event_scale();
    m_PtHat = genEventScale;
    
    int gcnt = 0; 
    m_nGen = 0;
    int nbs=0;
    
    for (HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();p != genEvent->particles_end(); ++p) {
      
      
      int mixed = -1;  // mixed is: -1 = unmixed
      GenVertex* endvert = (*p)->end_vertex(); 
      GenVertex* prodvert = (*p)->production_vertex();
      
      if (endvert && prodvert) {
	
	for ( GenVertex::particles_in_const_iterator p2 = prodvert->particles_in_const_begin(); p2 != prodvert->particles_in_const_end(); ++p2 ) { 
	  if ( (*p)->pdg_id() + (*p2)->pdg_id() == 0) {
	    mixed = 1;
	  }
	}
	
	for ( GenVertex::particles_out_const_iterator p22 = endvert->particles_out_const_begin(); p22 != endvert->particles_out_const_end(); ++p22 ) {
	  if ( (*p)->pdg_id() + (*p22)->pdg_id() == 0) mixed = 0;
	}
      }
      m_pdgId[gcnt]= (*p)->pdg_id();
      m_Number[gcnt]=(*p)->barcode() - 1;
      m_Status[gcnt]=(*p)->status();
      m_mcPx[gcnt]= (*p)->momentum().x();
      m_mcPy[gcnt]= (*p)->momentum().y();
      m_mcPz[gcnt]= (*p)->momentum().z();
      m_mcE[gcnt]= (*p)->momentum().e();
      m_mcPt[gcnt]= (*p)->momentum().perp();
      m_mcEta[gcnt]= (*p)->momentum().pseudoRapidity();
      m_mcPhi[gcnt]= (*p)->momentum().phi();
      
      GenVertex* pVertex = (*p)->end_vertex();
      if (0 != pVertex) {
	m_mcVx[gcnt] = pVertex->position().x(); 
	m_mcVy[gcnt] = pVertex->position().y();
	m_mcVz[gcnt] = pVertex->position().z(); 
      } else {
	m_mcVx[gcnt] = 9999.;
	m_mcVy[gcnt] = 9999.;
	m_mcVz[gcnt] = 9999.;
      }
      
      ++gcnt; 
      m_nGen = gcnt;
      TLorentzVector pmugenp;
      TLorentzVector pmugenm;
      TLorentzVector pKgenp;
      TLorentzVector pKgenm;
      int nBgen = 0;
      
      if ( abs((*p)->pdg_id()) == 531 && 0 != pVertex) {  // B_s       
	if (mixed!=0){
	  nbs++;
	}
	
	for ( GenVertex::particles_out_const_iterator ap = pVertex->particles_out_const_begin(); ap != pVertex->particles_out_const_end(); ++ap ) {
	  if ((*ap)->pdg_id() == 13) {  
	    pmugenp.SetPxPyPzE((*ap)->momentum().px(), (*ap)->momentum().py(),
			       (*ap)->momentum().pz(), (*ap)->momentum().e());
	  } else if ((*ap)->pdg_id() == -13) {  
	  pmugenm.SetPxPyPzE((*ap)->momentum().px(), (*ap)->momentum().py(),
			     (*ap)->momentum().pz(), (*ap)->momentum().e());
	  } else if ( (*ap)->pdg_id() == 333) {  
	    GenVertex* phiVertex = (*ap)->end_vertex();
	    for ( GenVertex::particles_out_const_iterator bp = phiVertex->particles_out_const_begin(); bp != phiVertex->particles_out_const_end(); ++bp ) {
	      if ( (*bp)->pdg_id() == 321) {  
		pKgenp.SetPxPyPzE((*bp)->momentum().px(), (*bp)->momentum().py(),
				  (*bp)->momentum().pz(), (*bp)->momentum().e());
	      } else if ( (*bp)->pdg_id() == -321) {  
		pKgenm.SetPxPyPzE((*bp)->momentum().px(), (*bp)->momentum().py(),
				  (*bp)->momentum().pz(), (*bp)->momentum().e());
	    }
	    }
	  }	
	}
	
	TLorentzVector pphigen = pKgenm + pKgenp;
	m_phiGenMass[nBgen] = pphigen.M();
	TLorentzVector pmumugen = pmugenp + pmugenm;
	m_mumuGenMass[nBgen] = pmumugen.M();
	TLorentzVector pbsgen = pphigen + pmumugen;
	m_bsGenMass[nBgen] = pbsgen.M();
	++nBgen; 
	m_nbGen = nBgen;
	
      }
    }

    m_nBs = nbs;

  } catch (...) {;}
}


void myAnalyzer::fillRecTracks(const edm::Event &iEvent) {
  // -- track collection
  edm::Handle<reco::TrackCollection> tracks;
  // edm::Handle<reco::CandidateCollection> tracks;
  iEvent.getByLabel(fTracksLabel.c_str(), tracks);  
  
  // -- track association
  m_nTrk = 0;
  int ntrk = 0;

  for (unsigned int i = 0; i < tracks->size(); ++i)
    { 
    reco::TrackRef rTrack(tracks, i);
    reco::Track track(*rTrack);
    
  // for (CandidateCollection::const_iterator tkCand = tracks->begin() ;
  //  tkCand !=tracks->end(); tkCand++ ){

    // const Track * track = dynamic_cast<const Track *> (&(*tkCand));

    // m_trkIndex[i]= rTrack.index();
    m_trkPt[ntrk]= track.pt(); 
    m_trkEta[ntrk]= track.eta();
    m_trkPhi[ntrk]= track.phi();
    m_trkIPt[ntrk]= track.d0();
    m_trkIPtErr[ntrk]= track.d0Error();
    m_trkIPz[ntrk]= track.dz();
    m_trkIPzErr[ntrk]= track.dzError();
    m_trkCharge[ntrk]= track.charge();
    m_trkChi2[ntrk]= track.chi2();
    m_trkDof[ntrk]= track.ndof();
    m_trkNHit[ntrk]= track.numberOfValidHits(); 

    ntrk++;
    m_nTrk = ntrk;
  }
}

void myAnalyzer::fillMuonBlock(const edm::Event &iEvent) {

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(fMuonLabel, muons);

  int munr = 0;
  m_nMu = 0;

  for (MuonCollection::const_iterator muonCand = muons->begin() ;
       muonCand !=muons->end(); muonCand++ ){

    if (fMuonLabel.label() == "muons") {

      // FullSim --> global muons
      m_muIndex[munr]= munr;
      m_muPt[munr]= muonCand->combinedMuon()->pt();
      m_muEta[munr]= muonCand->combinedMuon()->eta();
      m_muPhi[munr]= muonCand->combinedMuon()->phi();
      m_muIPt[munr]= muonCand->combinedMuon()->d0();
      m_muIPz[munr]= muonCand->combinedMuon()->dz();
      m_muChi2[munr]= muonCand->combinedMuon()->chi2();
      m_muCharge[munr]= muonCand->combinedMuon()->charge();
      m_muDof[munr]= muonCand->combinedMuon()->ndof();
      m_muNHit[munr]= muonCand->combinedMuon()->numberOfValidHits();

    } else {

      // FastSim --> parameterized tracks only
      m_muIndex[munr]= munr;
      m_muPt[munr]= muonCand->track()->pt();
      m_muEta[munr]= muonCand->track()->eta();
      m_muPhi[munr]= muonCand->track()->phi();
      m_muIPt[munr]= muonCand->track()->d0();
      m_muIPz[munr]= muonCand->track()->dz();
      m_muChi2[munr]= muonCand->track()->chi2();
      m_muCharge[munr]= muonCand->track()->charge();
      m_muDof[munr]= muonCand->track()->ndof();
      m_muNHit[munr]= muonCand->track()->numberOfValidHits();

    }
      
    munr++;
    m_nMu = munr;
  }
}

void myAnalyzer::fillJets(const edm::Event &iEvent) {

   edm::Handle<GenJetCollection> jetsgen;
   iEvent.getByLabel(fGenJetsLabel, jetsgen);

   int jtgnr = 0;
   m_nGJet = 0;

   for (unsigned int j = 0; j < jetsgen->size(); j++) {
     m_gjtEne[jtgnr] = (*jetsgen)[j].energy();
     m_gjtPhi[jtgnr] = (*jetsgen)[j].phi();
     m_gjtEta[jtgnr] = (*jetsgen)[j].eta();
     jtgnr++;
     m_nGJet = jtgnr;
   }

   edm::Handle<CaloJetCollection> jets;
   iEvent.getByLabel(fJetsLabel, jets);

   int jtnr = 0;
   m_nJet = 0;

   for (unsigned int i = 0; i < jets->size(); i++) {
     m_jtEne[jtnr] = (*jets)[i].energy();
     m_jtPhi[jtnr] = (*jets)[i].phi();
     m_jtEta[jtnr] = (*jets)[i].eta();
     jtnr++;
     m_nJet = jtnr;
   }
}

void myAnalyzer::fillTrigger(const edm::Event &iEvent) {
   edm::Handle<TriggerResults> trh;
   iEvent.getByLabel(fHLTLabel, trh);

   if (&trh) {
     int ntrigs = (*trh).size();
     if (ntrigs==0) {
       std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;
     }
     
     fTriggerNames.init( *trh );
    
     // Event 1 --> define root branches
     if (HltEvtCnt == 0) {
       for (int itrig = 0; itrig != ntrigs; ++itrig){
	 TString trigName = fTriggerNames.triggerName( itrig );
	 theTree->Branch(trigName, trigflag+itrig, trigName+"/I");
       }
       HltEvtCnt++;
     }

     // ...Fill the corresponding accepts in branch-variables
     for (int itrig = 0; itrig != ntrigs; ++itrig) {
       
       string trigName = fTriggerNames.triggerName( itrig );
       bool accept = (*trh).accept( itrig );
       
       if (accept) {
	 trigflag[itrig] = 1;
       } else {
	 trigflag[itrig] = 0;
       }
     }      
   }
}

void myAnalyzer::fillCands(const edm::Event &iEvent) {

  int phimumupairs = 0;
  m_bsCand = 0;
  m_bsVtx = 0;

  // get primary vertices  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(fVertexLabel,vertices);
  const Vertex * primaryVertex;
  if (vertices->size() > 0){
    primaryVertex = &(*vertices->begin());
    m_PVx = primaryVertex->x();
    m_PVy = primaryVertex->y();
    m_PVz = primaryVertex->z();
  } else {
    primaryVertex = 0;
    m_PVx = 9999.;
    m_PVy = 9999.;
    m_PVz = 9999.;
  }

  edm::Handle<reco::CandidateCollection> bCandCollection;
  iEvent.getByLabel(fBCandLabel,bCandCollection);
  // edm::Handle<reco::VertexCompositeCandidateCollection> bVtxCollection;
  // iEvent.getByLabel(fBVtxLabel,bVtxCollection);
 
  cout << "Found: " << bCandCollection->size() << " reconstructed Bs" << "\n";
  // cout << "Found: " << bVtxCollection->size() << " vertexed Bs" << "\n";

  // TEST
  edm::Handle<reco::CandidateCollection> allTrackCollection;
  iEvent.getByLabel("allTracks",allTrackCollection);
  cout << "Found: " << allTrackCollection->size() << " good tracks" << "\n";
  edm::Handle<reco::CandidateCollection> allMuonCollection;
  iEvent.getByLabel("allMuons",allMuonCollection);
  cout << "Found: " << allMuonCollection->size() << " good muons" << "\n";
  edm::Handle<reco::CandidateCollection> mumuCollection;
  iEvent.getByLabel("CandMuMu",mumuCollection);
  cout << "Found: " << mumuCollection->size() << " mu mu pairs" << "\n";
  edm::Handle<reco::CandidateCollection> phiCollection;
  iEvent.getByLabel("CandPhi",phiCollection);
  cout << "Found: " << phiCollection->size() << " phi candidates" << "\n"; 
  //

  float theMaxVar = 0.;
  if (storeTheBest && whichBest == "PhiMass") theMaxVar = 100000.;
  for (CandidateCollection::const_iterator bit=(*bCandCollection).begin();bit!=(*bCandCollection).end();bit++) {
				   
    if ((*bit).numberOfDaughters() != 2){
      continue;
    }

    // get any mu-mu candidate
    const Candidate *mumu = bit->daughter(0);  

    TransientTrack muPlusTk;				
    TransientTrack muMinusTk;
    for (Candidate::const_iterator muit = mumu->begin(); muit!=mumu->end(); muit++) {	
      const Muon *mucand = dynamic_cast<const Muon *> (&(*muit));
      TrackRef mutk = mucand->track();
      TransientTrack muttkp   = (*theB).build(&mutk);			   
      if (muttkp.charge() > 0) {
	muPlusTk = muttkp; 
      } else {
        muMinusTk =  muttkp;
      }
    } 
   
    // get any phi candidate
    const Candidate *thephi = bit->daughter(1);
    TransientTrack kPlusTk;				
    TransientTrack kMinusTk;		    				
    for (Candidate::const_iterator kit = thephi->begin(); kit!=thephi->end(); kit++){
      const RecoCandidate* kcand = dynamic_cast<const RecoCandidate*> (&(*kit));
      TrackRef ktk = kcand->track();
      TransientTrack kttkp   = (*theB).build(&ktk);
      if (kttkp.charge() > 0 ){
	kPlusTk = kttkp;
      } else {
        kMinusTk =  kttkp;
      }
    } 

    float tmpbsmass = (*bit).mass();
    float tmpmumumass = mumu->mass();
    float tmpphimass = thephi->mass();
    float tmpisolationdr12 = this->isolation(iEvent, &(*bit), 1.2);
    float tmpisolationdr1 = this->isolation(iEvent, &(*bit), 1.0);
    float tmpisolationdr07 = this->isolation(iEvent, &(*bit), 0.7);
    float tmpisolationdr05 = this->isolation(iEvent, &(*bit), 0.5);

    // VERTICE (lotte)
    KinematicParticleFactoryFromTransientTrack pFactory;
    
    // The mass of a muon and the insignificant mass sigma to avoid singularities in the covariance matrix.
    ParticleMass muon_mass = 0.1056583;
    ParticleMass kaon_mass = 0.493677;
    float muon_sigma = 0.0000000001;
    float kaon_sigma = 0.000016;
			
    float chi = 0.;
    float ndf = 0.;
    
    // making particles
    vector<RefCountedKinematicParticle> allParticles;
    allParticles.push_back(pFactory.particle (kPlusTk,kaon_mass,chi,ndf,kaon_sigma));
    allParticles.push_back(pFactory.particle (kMinusTk,kaon_mass,chi,ndf,kaon_sigma));
    allParticles.push_back(pFactory.particle (muPlusTk,muon_mass,chi,ndf,muon_sigma));
    allParticles.push_back(pFactory.particle (muMinusTk,muon_mass,chi,ndf,muon_sigma));
    
    // creating the constraint for the phi mass
    ParticleMass m_phi = 1.01946;
    MultiTrackKinematicConstraint *phi_c = new TwoTrackMassKinematicConstraint(m_phi);
    
    // fit to the vertex
    KinematicConstrainedVertexFitter kcvFitter;
    RefCountedKinematicTree myTree = kcvFitter.fit(allParticles, phi_c);
    cout << "Global vertex fit done\n";
    
    myTree->movePointerToTheTop();
    RefCountedKinematicParticle newBs = myTree->currentParticle();
    RefCountedKinematicVertex bVertex = myTree->currentDecayVertex();

    // get fit properties
    float fchi2 = newBs->chiSquared();
    int fndof =(int)newBs->degreesOfFreedom();
    float tmpprob = TMath::Prob(fchi2, fndof);
    
    // get fitted candidates
    float tmpbsmassvtx = newBs->currentState().mass(); 
    float tmpbsptvtx = newBs->currentState().globalMomentum().perp();
    float tmpbsetavtx = newBs->currentState().globalMomentum().eta();
    float tmpbsphivtx = newBs->currentState().globalMomentum().phi();

    float tmpphimassvtx = -999.;
    float tmpmumumassvtx = -999.;
    float tmpmumuDR = -999.;

    vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
    if(bs_children.size() == 4) {

       TLorentzVector pmufitp;
       TLorentzVector pmufitm;
       TLorentzVector pKfitp;
       TLorentzVector pKfitm;
       pKfitp.SetPxPyPzE(bs_children[0]->currentState().globalMomentum().x(),
			 bs_children[0]->currentState().globalMomentum().y(),
                         bs_children[0]->currentState().globalMomentum().z(),
			 sqrt( pow (bs_children[0]->currentState().globalMomentum().mag(),2) + pow (bs_children[0]->currentState().mass(),2) ));
       pKfitm.SetPxPyPzE(bs_children[1]->currentState().globalMomentum().x(),
			 bs_children[1]->currentState().globalMomentum().y(),
                         bs_children[1]->currentState().globalMomentum().z(),
			 sqrt( pow (bs_children[1]->currentState().globalMomentum().mag(),2) + pow (bs_children[1]->currentState().mass(),2) ));
       pmufitp.SetPxPyPzE(bs_children[2]->currentState().globalMomentum().x(),
			 bs_children[2]->currentState().globalMomentum().y(),
                         bs_children[2]->currentState().globalMomentum().z(),
			 sqrt( pow (bs_children[2]->currentState().globalMomentum().mag(),2) + pow (bs_children[2]->currentState().mass(),2) ));
       pmufitm.SetPxPyPzE(bs_children[3]->currentState().globalMomentum().x(),
			 bs_children[3]->currentState().globalMomentum().y(),
                         bs_children[3]->currentState().globalMomentum().z(),
			 sqrt( pow (bs_children[3]->currentState().globalMomentum().mag(),2) + pow (bs_children[3]->currentState().mass(),2) ));
       
       TLorentzVector pphifit = pKfitm + pKfitp;
       tmpphimassvtx = pphifit.M();
       TLorentzVector pmumufit = pmufitp + pmufitm;
       tmpmumumassvtx = pmumufit.M();

       float deta = bs_children[3]->currentState().globalMomentum().eta() - bs_children[2]->currentState().globalMomentum().eta();
       float dphi = bs_children[3]->currentState().globalMomentum().phi() - bs_children[2]->currentState().globalMomentum().phi();
       if (dphi > pai) dphi = 2*pai - dphi;
       tmpmumuDR = sqrt(deta*deta + dphi*dphi);

    } else { cout << "B_s had " << bs_children.size() << " children" << endl; }

    // get decay lengths and pointing
    float tmpsvx = bVertex->vertexState().position().x();
    float tmpsvy = bVertex->vertexState().position().y();
    float tmpsvz = bVertex->vertexState().position().z();

    float tmpdist3d = -999.;
    float tmpdist3derr = -999.;
    float tmpdist2d = -999.;
    float tmpdist2derr = -999.;
    float tmppointing = -999.;
    float tmpcospointing = -999.;

    if (primaryVertex) {
      VertexDistance3D vdist3d;
      Measurement1D dist3d = vdist3d.distance(bVertex->vertexState(), *primaryVertex);
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance(bVertex->vertexState(), *primaryVertex);
      tmpdist3d = dist3d.value();
      tmpdist3derr = dist3d.error();
      tmpdist2d = distXY.value();
      tmpdist2derr = distXY.error();
      
      TVector3 vecdist;
      TVector3 vecmom;
      vecdist.SetXYZ(tmpsvx - m_PVx,tmpsvy - m_PVy,tmpsvz - m_PVz);
      vecmom.SetXYZ(newBs->currentState().globalMomentum().x(),            
		    newBs->currentState().globalMomentum().y(),
		    newBs->currentState().globalMomentum().z());
      tmppointing = vecdist.Angle( vecmom );
      tmpcospointing = cos(tmppointing);
    }

    if (!storeTheBest) {
      
      m_bsMass[phimumupairs] = tmpbsmass;
      m_mumuMass[phimumupairs] = tmpmumumass;
      m_phiMass[phimumupairs] = tmpphimass;
      m_bsIsoDR12[phimumupairs] = tmpisolationdr12;
      m_bsIsoDR1[phimumupairs] = tmpisolationdr1;
      m_bsIsoDR07[phimumupairs] = tmpisolationdr07;
      m_bsIsoDR05[phimumupairs] = tmpisolationdr05;
      m_bsMassVtx[phimumupairs] = tmpbsmassvtx;
      m_bsPtVtx[phimumupairs] = tmpbsptvtx;
      m_bsEtaVtx[phimumupairs] = tmpbsetavtx;
      m_bsPhiVtx[phimumupairs] = tmpbsphivtx;
      m_mumuMassVtx[phimumupairs] = tmpmumumassvtx;
      m_mumuDRVtx[phimumupairs] = tmpmumuDR;
      m_phiMassVtx[phimumupairs] = tmpphimassvtx;
      m_chi2Vtx[phimumupairs] = fchi2;
      m_probVtx[phimumupairs] = tmpprob;
      m_SVx[phimumupairs] = tmpsvx;
      m_SVy[phimumupairs] = tmpsvy;
      m_SVz[phimumupairs] = tmpsvz;
      m_decLen2D[phimumupairs] = tmpdist2d;
      m_decSigma2D[phimumupairs]= tmpdist2derr;
      m_decLen3D[phimumupairs] = tmpdist3d;
      m_decSigma3D[phimumupairs] = tmpdist3derr;
      m_bsPointing[phimumupairs] = tmppointing;
      m_bsCosPointing[phimumupairs] = tmpcospointing;
      ++phimumupairs;
      m_bsCand = phimumupairs;
      m_bsVtx = phimumupairs;

    } else {
      
      bool compare = false;
      float thisVar;
      if (whichBest == "KsWithHighestPt") {

	const Candidate *theKtrackp = thephi->daughter(0);
	const Candidate *theKtrackm = thephi->daughter(1);
	thisVar = theKtrackp->pt() + theKtrackm->pt();
        compare = (thisVar > theMaxVar);

      } else if (whichBest == "PhiMass") {
       
        thisVar = fabs(thephi->mass() - 1.01946);
	compare = (thisVar < theMaxVar);

      }

      if (compare) {
	m_bsMass[0] = tmpbsmass;
	m_mumuMass[0] = tmpmumumass;
	m_phiMass[0] = tmpphimass;
        m_bsIsoDR12[0] = tmpisolationdr12;
        m_bsIsoDR1[0] = tmpisolationdr1;
        m_bsIsoDR07[0] = tmpisolationdr07;
        m_bsIsoDR05[0] = tmpisolationdr05;
        m_bsMassVtx[0] = tmpbsmassvtx;
	m_bsPtVtx[0] = tmpbsptvtx;
	m_bsEtaVtx[0] = tmpbsetavtx;
	m_bsPhiVtx[0] = tmpbsphivtx;
	m_mumuDRVtx[0] = tmpmumuDR;
	m_mumuMassVtx[0] = tmpmumumassvtx;
	m_phiMassVtx[0] = tmpphimassvtx;
	m_chi2Vtx[0] = fchi2;
	m_probVtx[0] = tmpprob;
        m_SVx[0] = tmpsvx;
	m_SVy[0] = tmpsvy;
	m_SVz[0] = tmpsvz;
	m_decLen2D[0] = tmpdist2d;
	m_decSigma2D[0]= tmpdist2derr;
	m_decLen3D[0] = tmpdist3d;
	m_decSigma3D[0] = tmpdist3derr;
	m_bsPointing[0] = tmppointing;
	m_bsCosPointing[0] = tmpcospointing;
	m_bsCand = 1;
        m_bsVtx = 1;
	thisVar = theMaxVar;
      }
    }     
  }

}

float  myAnalyzer::isolation(const edm::Event &iEvent, const Candidate* aCand, double cone) {

  float pai = 3.1415926;

  edm::Handle<reco::TrackCollection> thetracks;
  iEvent.getByLabel(fTracksLabel.c_str(), thetracks);  
  
  const Candidate *oneDaug = aCand->daughter(0);   
  const Candidate *anotherDaug = aCand->daughter(1);
  const Candidate *gDaug[4];
  gDaug[0] = oneDaug->daughter(0);
  gDaug[1] = oneDaug->daughter(1); 
  gDaug[2] = anotherDaug->daughter(0);
  gDaug[3] = anotherDaug->daughter(1);

  float thisIso = 0.;
  for (unsigned int i = 0; i < thetracks->size(); ++i){

    bool itOverlaps = false; 
    reco::TrackRef rTrack(thetracks, i);
    reco::Track track(*rTrack);

    float deltaeta = track.eta()-aCand->eta();
    float deltaphi = fabs(track.phi()-aCand->phi());
    if (deltaphi > pai) deltaphi = 2*pai - deltaphi;
    float dr = sqrt( pow(deltaphi,2) + pow(deltaeta,2) );

    for (unsigned int igdaug = 0; igdaug < 4; igdaug++) {
      if ( fabs(track.phi() - gDaug[igdaug]->phi()) < 0.0001) itOverlaps = true ;
      // cout << "TRY:" << endl;
      // cout << "Track phi = " << track.phi() << " Daug phi = " << gDaug[igdaug]->phi() << " / Track eta = " << track.eta() << " Daug phi = " << gDaug[igdaug]->eta() << endl;
    }
    if (!itOverlaps && dr < cone) thisIso += track.pt();
  }
  thisIso = aCand->pt() / (aCand->pt() + thisIso);
  return thisIso;
}

void  myAnalyzer::beginJob(const EventSetup& setup) {
  // get the builder to transform tracks to TransientTrack
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);	
}

// ------------ method called once each job just after ending the event loop  ------------
void  myAnalyzer::endJob() {
  theFile->cd();
  theTree->Write();
  theFile->Close();
  delete theFile;
}

//define this as a plug-in
//DEFINE_FWK_MODULE(myAnalyzer);
