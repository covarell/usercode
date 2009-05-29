#include <TStyle.h>
#include <TCanvas.h>
#include "TLine.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string.h>

#include <math.h>
#include <map>
#include <vector>

#include "finalJPsiAnalysis.hh"

using namespace std;

finalJPsiAnalysis::finalJPsiAnalysis(TTree *tree) 
  : JPsiTreeBase(tree) {
}

finalJPsiAnalysis::~finalJPsiAnalysis(){ } 

void finalJPsiAnalysis::Loop() {

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // booking histos
  bookHistos();

  // loop over events
  std::cout << "Number of entries = " << nentries << std::endl;

  // counters
  int totalEvents      = 0;
                
  for (int jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    totalEvents++;

    // cout << "Passed? " << (int)HLTBits_accept[HLTbits[0]] << endl;
    // if (!(HLTBits_accept[HLTbits[0]])) continue;

    bool aRealJpsiEvent = false;
    for (int imugen=0; imugen<Mc_mu_size; imugen++) {
      if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
    }
    if (aRealJpsiEvent) continue;

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {
      TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
      float theMass = theQQ4mom->M();
      if (Reco_QQ_sign[iqq] == 0) {
	if (HLTBits_accept[0]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_passmu3->Fill(theMass);
          if (Reco_QQ_type[iqq] == 0 && theQQ4mom->Perp() < 6.0) QQMass2GlobPT6_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1 && theQQ4mom->Perp() < 6.0) QQMass1Glob1TrkPT6_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3 && theQQ4mom->Perp() < 6.0) QQMass1Glob1CalPT6_passmu3->Fill(theMass);
	}
	if (HLTBits_accept[1]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_passmu5->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_passmu5->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_passmu5->Fill(theMass);
          if (Reco_QQ_type[iqq] == 0 && theQQ4mom->Perp() < 6.0) QQMass2GlobPT6_passmu5->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1 && theQQ4mom->Perp() < 6.0) QQMass1Glob1TrkPT6_passmu5->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3 && theQQ4mom->Perp() < 6.0) QQMass1Glob1CalPT6_passmu5->Fill(theMass); 
	}
	if (HLTBits_accept[2]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_passmu9->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_passmu9->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_passmu9->Fill(theMass);
	}
	if (HLTBits_accept[4]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_pass2isomu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_pass2isomu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_pass2isomu3->Fill(theMass);
	} 
	if (HLTBits_accept[5]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_pass2mu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_pass2mu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_pass2mu3->Fill(theMass);
	}
	if (HLTBits_accept[6]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_pass2mu3jpsi->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_pass2mu3jpsi->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_pass2mu3jpsi->Fill(theMass);
	}
	if (HLTBits_accept[7]) {
	  if (Reco_QQ_type[iqq] == 0) QQMass2Glob_pass2mu3ups->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) QQMass1Glob1Trk_pass2mu3ups->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) QQMass1Glob1Cal_pass2mu3ups->Fill(theMass);
	}
      } else {
	if (HLTBits_accept[0]) {
	  if (Reco_QQ_type[iqq] == 0) WSMass2Glob_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) WSMass1Glob1Trk_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) WSMass1Glob1Cal_passmu3->Fill(theMass);
	}
        if (HLTBits_accept[5]) {
	  if (Reco_QQ_type[iqq] == 0) WSMass2Glob_pass2mu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) WSMass1Glob1Trk_pass2mu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) WSMass1Glob1Cal_pass2mu3->Fill(theMass);
	}
      }
    }
  }

  saveHistos();

  drawPlots();

  // summary
  cout << "total number of events = "              << totalEvents     << endl;

} // end of program


void finalJPsiAnalysis::bookHistos() {

  QQMass2Glob_passmu3              = new TH1F("QQMass2Glob_passmu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_passmu3          = new TH1F("QQMass1Glob1Trk_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_passmu3          = new TH1F("QQMass1Glob1Cal_passmu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2GlobPT6_passmu3              = new TH1F("QQMass2GlobPT6_passmu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1TrkPT6_passmu3          = new TH1F("QQMass1Glob1TrkPT6_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1CalPT6_passmu3          = new TH1F("QQMass1Glob1CalPT6_passmu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  WSMass2Glob_passmu3              = new TH1F("WSMass2Glob_passmu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  WSMass1Glob1Trk_passmu3          = new TH1F("WSMass1Glob1Trk_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  WSMass1Glob1Cal_passmu3          = new TH1F("WSMass1Glob1Cal_passmu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2Glob_passmu5              = new TH1F("QQMass2Glob_passmu5",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_passmu5          = new TH1F("QQMass1Glob1Trk_passmu5",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_passmu5          = new TH1F("QQMass1Glob1Cal_passmu5",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.); 
  QQMass2Glob_passmu9              = new TH1F("QQMass2Glob_passmu9",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_passmu9          = new TH1F("QQMass1Glob1Trk_passmu9",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_passmu9          = new TH1F("QQMass1Glob1Cal_passmu9",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2GlobPT6_passmu5              = new TH1F("QQMass2GlobPT6_passmu5",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1TrkPT6_passmu5          = new TH1F("QQMass1Glob1TrkPT6_passmu5",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1CalPT6_passmu5          = new TH1F("QQMass1Glob1CalPT6_passmu5",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2Glob_pass2mu3              = new TH1F("QQMass2Glob_pass2mu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_pass2mu3          = new TH1F("QQMass1Glob1Trk_pass2mu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_pass2mu3          = new TH1F("QQMass1Glob1Cal_pass2mu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  WSMass2Glob_pass2mu3              = new TH1F("WSMass2Glob_pass2mu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  WSMass1Glob1Trk_pass2mu3          = new TH1F("WSMass1Glob1Trk_pass2mu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  WSMass1Glob1Cal_pass2mu3          = new TH1F("WSMass1Glob1Cal_pass2mu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2Glob_pass2isomu3              = new TH1F("QQMass2Glob_pass2isomu3",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_pass2isomu3          = new TH1F("QQMass1Glob1Trk_pass2isomu3",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_pass2isomu3          = new TH1F("QQMass1Glob1Cal_pass2isomu3",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2Glob_pass2mu3jpsi              = new TH1F("QQMass2Glob_pass2mu3jpsi",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_pass2mu3jpsi          = new TH1F("QQMass1Glob1Trk_pass2mu3jpsi",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_pass2mu3jpsi          = new TH1F("QQMass1Glob1Cal_pass2mu3jpsi",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
  QQMass2Glob_pass2mu3ups              = new TH1F("QQMass2Glob_pass2mu3ups",  "Invariant mass (2 global muons)", 100, 0.,15.);
  QQMass1Glob1Trk_pass2mu3ups          = new TH1F("QQMass1Glob1Trk_pass2mu3ups",  "Invariant mass (1 global + 1 tracker muon)", 100, 0.,15.);
  QQMass1Glob1Cal_pass2mu3ups          = new TH1F("QQMass1Glob1Cal_pass2mu3ups",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
}

void finalJPsiAnalysis::drawPlots() {

  // TCanvas c; c.cd();
  // EleHistoGt4_size              -> Draw();  // c.Print("EleGt4.eps");
  
}

void finalJPsiAnalysis::saveHistos() {

  TFile fOut("Outfile_histo.root", "RECREATE");
  fOut.cd();
  
  QQMass2Glob_passmu3        -> Write();        
  QQMass1Glob1Trk_passmu3    -> Write(); 
  QQMass1Glob1Cal_passmu3    -> Write();
  QQMass2GlobPT6_passmu3        -> Write();        
  QQMass1Glob1TrkPT6_passmu3    -> Write(); 
  QQMass1Glob1CalPT6_passmu3    -> Write();
  WSMass2Glob_passmu3        -> Write();        
  WSMass1Glob1Trk_passmu3    -> Write(); 
  WSMass1Glob1Cal_passmu3    -> Write();
  QQMass2Glob_passmu5        -> Write();
  QQMass1Glob1Trk_passmu5    -> Write();
  QQMass1Glob1Cal_passmu5    -> Write();
  QQMass2GlobPT6_passmu5        -> Write();
  QQMass1Glob1TrkPT6_passmu5    -> Write();
  QQMass1Glob1CalPT6_passmu5    -> Write();
  QQMass2Glob_passmu9        -> Write();
  QQMass1Glob1Trk_passmu9    -> Write();
  QQMass1Glob1Cal_passmu9    -> Write();
  QQMass2Glob_pass2mu3       -> Write();
  QQMass1Glob1Trk_pass2mu3   -> Write();
  QQMass1Glob1Cal_pass2mu3   -> Write();
  WSMass2Glob_pass2mu3       -> Write();
  WSMass1Glob1Trk_pass2mu3   -> Write();
  WSMass1Glob1Cal_pass2mu3   -> Write();
  QQMass2Glob_pass2isomu3      -> Write();  
  QQMass1Glob1Trk_pass2isomu3  -> Write();  
  QQMass1Glob1Cal_pass2isomu3  -> Write();  
  QQMass2Glob_pass2mu3jpsi     -> Write();  
  QQMass1Glob1Trk_pass2mu3jpsi -> Write();  
  QQMass1Glob1Cal_pass2mu3jpsi -> Write();  
  QQMass2Glob_pass2mu3ups      -> Write();  
  QQMass1Glob1Trk_pass2mu3ups  -> Write();  
  QQMass1Glob1Cal_pass2mu3ups  -> Write();  
}			      
			       
			       
