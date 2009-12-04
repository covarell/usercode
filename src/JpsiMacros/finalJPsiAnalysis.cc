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

  /// SELECTION CUTS ///
 
  MIN_nhits_trk = 12;
  MAX_normchi2_trk = 5.0;
  MAX_normchi2_glb = 20.0;
  MIN_nhits_pixel = 2;
  MAX_d0_trk = 5.0;
  MAX_dz_trk = 20.0;
  MIN_vtxprob_jpsi = 0.001;

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
    if (jentry%50000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    totalEvents++;

    int theMCMatchedGlbMu1 = -1;
    int theMCMatchedGlbMu2 = -1;
    int theMCMatchedTrkMu = -1;
    int theMCMatchedCalMu = -1;

    // cout << "Passed? " << (int)HLTBits_accept[HLTbits[0]] << endl;

    // *** CUT ON TRIGGER ***
    // if (!(HLTBits_accept[0])) continue;

    passingHLTMu3->Fill(int(HLTBits_accept[0]));  
    trackSize->Fill(Reco_track_size);      
    gammaSize->Fill(Reco_gamma_size);      
    globalMuonSize->Fill(Reco_mu_glb_size); 
    trackerMuonSize->Fill(Reco_mu_trk_size);
    caloMuonSize->Fill(Reco_mu_cal_size);  

    bool aRealJpsiEvent = false;
    unsigned int nGlobMuMatched = 0;

    //
    // Find MC matched muons (all QQ categories)
    //
    for (int imugen=0; imugen<Mc_mu_size; imugen++) {
      // cout << "Mc_mu_size = " << Mc_mu_size << endl;
      if (Mc_mumoth_id[imugen] == 443) {
        TLorentzVector *theMcMumom = (TLorentzVector*)Mc_mu_4mom->At(imugen);
        cout << "Found a MC muon " << Mc_mu_id[imugen] << " " << theMcMumom->Perp() << " " << theMcMumom->Eta() << endl;
	aRealJpsiEvent = true;
        for (int imugl=0; imugl<Reco_mu_glb_size; imugl++) {
          // cout << "Reco_mu_glb_size = " << Reco_mu_glb_size << endl;
          TLorentzVector *theGlMumom = (TLorentzVector*)Reco_mu_glb_4mom->At(imugl);
          hMcRecoGlobMuDeltaR->Fill(deltaR(theMcMumom,theGlMumom));
          if (deltaR(theMcMumom,theGlMumom) < 0.03) {
            if (nGlobMuMatched == 0) {
	      theMCMatchedGlbMu1 = imugl;
              nGlobMuMatched++;
	    } else {
	      theMCMatchedGlbMu2 = imugl;
	      break;
	    }
	  }
	}
        for (int imutr=0; imutr<Reco_mu_trk_size; imutr++) {
          // cout << "Reco_mu_trk_size = " << Reco_mu_trk_size << endl;
          TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(imutr);
          hMcRecoTrkMuDeltaR->Fill(deltaR(theMcMumom,theTrMumom));
          if (deltaR(theMcMumom,theTrMumom) < 0.03) {
	    theMCMatchedTrkMu = imutr;      break;
	  }
	}
         for (int imuca=0; imuca<Reco_mu_cal_size; imuca++) {
	   // cout << "Reco_mu_cal_size = " << Reco_mu_cal_size << endl;
          TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(imuca);
          hMcRecoCalMuDeltaR->Fill(deltaR(theMcMumom,theCaMumom));
          if (deltaR(theMcMumom,theCaMumom) < 0.03) {
	    theMCMatchedCalMu = imuca;       break;
	  }
	}
      }      
    }
    /* cout << "theMCMatchedGlbMu1 = " << theMCMatchedGlbMu1 << endl;
    cout << "theMCMatchedGlbMu2 = " << theMCMatchedGlbMu2 << endl;
    cout << "theMCMatchedTrkMu = " << theMCMatchedTrkMu << endl;
    cout << "theMCMatchedCalMu = " << theMCMatchedCalMu << endl;
    cout << "************************************** " << endl; */

    genQQSize->Fill(Mc_QQ_size);
    if (Mc_QQ_size > 0) cout << "Found a MC J/psi " << Mc_QQmoth_id[0] << endl;
    if (aRealJpsiEvent) {
      TLorentzVector *theQQMcmom = (TLorentzVector*)Mc_QQ_4mom->At(0);
      genQQPt->Fill(theQQMcmom->Perp());
    }

    //
    // MC matched global mu (before cuts)
    //
    for (int imugl=0; imugl<Reco_mu_glb_size; imugl++) {
      TLorentzVector *theGlMumom = (TLorentzVector*)Reco_mu_glb_4mom->At(imugl);
      globalMuonPt->Fill(theGlMumom->Perp());
      if (theMCMatchedGlbMu1 == imugl || theMCMatchedGlbMu2 == imugl) {
        hMcRightAllMuIso->Fill(Reco_mu_glb_iso[imugl]);
	hMcRightGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
        hMcRightGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
        hMcRightGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
        if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	  hMcRightGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);
      } else {
        hMcWrongGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
        hMcWrongGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
        hMcWrongGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
        if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	  hMcWrongGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);

      }
    }

    //
    // MC matched calo mu (before cuts)
    //
    for (int imuca=0; imuca<Reco_mu_cal_size; imuca++) {
      
      TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(imuca);
      caloMuonPt->Fill(theCaMumom->Perp());
      if (theMCMatchedCalMu == imuca) {
	hMcRightCalMuPt->Fill(theCaMumom->Perp());
        hMcRightCalMuNhits->Fill(Reco_mu_cal_nhitstrack[imuca]);
        hMcRightCalMuChi2->Fill(Reco_mu_cal_normChi2[imuca]);
	hMcRightCalMuCaloComp->Fill(Reco_mu_cal_caloComp[imuca]);
      } else {
	hMcWrongCalMuPt->Fill(theCaMumom->Perp());
        hMcWrongCalMuNhits->Fill(Reco_mu_cal_nhitstrack[imuca]);
        hMcWrongCalMuChi2->Fill(Reco_mu_cal_normChi2[imuca]);
	hMcWrongCalMuCaloComp->Fill(Reco_mu_cal_caloComp[imuca]);
      }
    }

    //
    // MC matched tracker mu (before cuts)
    //
    for (int imutr=0; imutr<Reco_mu_trk_size; imutr++) {
      //cout << "imutr = " << imutr << " TrkSize = " << Reco_mu_trk_size << endl;
      TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(imutr);
      trackerMuonPt->Fill(theTrMumom->Perp());
      if (theMCMatchedTrkMu == imutr) {
	hMcRightAllMuIso->Fill(Reco_mu_trk_iso[imutr]);
	hMcRightTrkMuPt->Fill(theTrMumom->Perp());
        hMcRightTrkMuNhits->Fill(Reco_mu_trk_nhitstrack[imutr]);
        hMcRightTrkMuChi2->Fill(Reco_mu_trk_normChi2[imutr]);
        hMcRightTrkMud0->Fill(Reco_mu_trk_d0[imutr]);
        hMcRightTrkMudz->Fill(Reco_mu_trk_dz[imutr]);
	hMcRightTrkBit4->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,4))/(int)pow(2,4));
	hMcRightTrkBit5->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,5))/(int)pow(2,5));
	hMcRightTrkBit8->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,8))/(int)pow(2,8));
	hMcRightTrkBit9->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,9))/(int)pow(2,9));
      } else {
	hMcWrongTrkMuPt->Fill(theTrMumom->Perp());
	hMcWrongTrkMuNhits->Fill(Reco_mu_trk_nhitstrack[imutr]);
        hMcWrongTrkMuChi2->Fill(Reco_mu_trk_normChi2[imutr]);
	hMcWrongTrkMud0->Fill(Reco_mu_trk_d0[imutr]);
        hMcWrongTrkMudz->Fill(Reco_mu_trk_dz[imutr]);
	hMcWrongTrkBit4->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,4))/(int)pow(2,4));
	hMcWrongTrkBit5->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,5))/(int)pow(2,5));
	hMcWrongTrkBit8->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,8))/(int)pow(2,8));  
	hMcWrongTrkBit9->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,9))/(int)pow(2,9));
      }
    }
    
  

    int bestQQ = theBestQQ();

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {
      TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
      float theMass = theQQ4mom->M();
      if (Reco_QQ_sign[iqq] == 0) {
	if (Reco_QQ_type[iqq] == 0) {
	  QQMass2Glob_passmu3->Fill(theMass);
          QQPt2Glob_passmu3->Fill(theQQ4mom->Perp());
          if ( (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu1) || (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu2) ) {
	      hMcRightGlbGlbMuMass->Fill(theMass);  
	    } else {
	      hMcWrongGlbGlbMuMass->Fill(theMass);
	    }
	  if (iqq == bestQQ) {
	    QQMass2Glob_best->Fill(theMass);
	    
	    //
	    // MC matched global mu (after cuts)
	    //
	    /* if ( (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu1) || (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu2) ) {
	       hMcRightGlbGlbMuVtxProb->Fill(Reco_QQ_probChi2[iqq]);  
	       } else {
	       hMcWrongGlbGlbMuVtxProb->Fill(Reco_QQ_probChi2[iqq]);
	       }
	       
	       int imugl = Reco_QQ_muhpt[iqq];
	       if (imugl == theMCMatchedGlbMu1 || imugl == theMCMatchedGlbMu2) {
	       hMcRightGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
	       hMcRightGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
	       hMcRightGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
	       if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	       hMcRightGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);
	       } else {
	       hMcWrongGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
	       hMcWrongGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
	       hMcWrongGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
	       if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	       hMcWrongGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);
	       }
	       
	       imugl = Reco_QQ_mulpt[iqq];
	       if (Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu2) {
	       hMcRightGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
	       hMcRightGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
	       hMcRightGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
	       if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	       hMcRightGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);
	       } else {
	       hMcWrongGlbMunPixHits->Fill(Reco_mu_glb_nhitsPixB[imugl] + Reco_mu_glb_nhitsPixE[imugl]);
	       hMcWrongGlbMud0->Fill(Reco_mu_glb_d0[imugl]);
	       hMcWrongGlbMudz->Fill(Reco_mu_glb_dz[imugl]);
	       if (Reco_mu_glb_nhitsPix1HitBE[imugl] > 0) 
	       hMcWrongGlbMuFirstLayer->Fill(Reco_mu_glb_nhitsPix1Hit[imugl]);
	       } */
	  }
	}
	
	if (Reco_QQ_type[iqq] == 1) {
	  QQMass1Glob1Trk_passmu3->Fill(theMass);
          QQPt1Glob1Trk_passmu3->Fill(theQQ4mom->Perp());
          if ( Reco_QQ_mulpt[iqq] == theMCMatchedTrkMu && (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2) ) {
	    hMcRightGlbTrkMuMass->Fill(theMass);  
	  } else {
	    hMcWrongGlbTrkMuMass->Fill(theMass);
	  }
	  if (iqq == bestQQ) {
	    QQMass1Glob1Trk_best->Fill(theMass);

	    //
	    // MC matched trk mu (after cuts)
	    //
	    /* if ( Reco_QQ_mulpt[iqq] == theMCMatchedTrkMu && (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2) ) {
	       hMcRightGlbTrkMuVtxProb->Fill(Reco_QQ_probChi2[iqq]);  
	       } else {
	       hMcWrongGlbTrkMuVtxProb->Fill(Reco_QQ_probChi2[iqq]);
	       }
	       
	       int imutr = Reco_QQ_mulpt[iqq];
	       TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(imutr);
	       if (theMCMatchedTrkMu == imutr) {
	       hMcRightTrkMuPt->Fill(theTrMumom->Perp());
	       hMcRightTrkMuNhits->Fill(Reco_mu_trk_nhitstrack[imutr]);
	       hMcRightTrkMuChi2->Fill(Reco_mu_trk_normChi2[imutr]);
	       hMcRightTrkMud0->Fill(Reco_mu_trk_d0[imutr]);
	       hMcRightTrkMudz->Fill(Reco_mu_trk_dz[imutr]);
	       hMcRightTrkBit4->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,4))/(int)pow(2,4));
	       hMcRightTrkBit5->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,5))/(int)pow(2,5));
	       hMcRightTrkBit8->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,8))/(int)pow(2,8));
	       hMcRightTrkBit9->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,9))/(int)pow(2,9));
	       } else {
	       hMcWrongTrkMuPt->Fill(theTrMumom->Perp());
	       hMcWrongTrkMuNhits->Fill(Reco_mu_trk_nhitstrack[imutr]);
	       hMcWrongTrkMuChi2->Fill(Reco_mu_trk_normChi2[imutr]);
	       hMcWrongTrkMud0->Fill(Reco_mu_trk_d0[imutr]);
	       hMcWrongTrkMudz->Fill(Reco_mu_trk_dz[imutr]);
	       hMcWrongTrkBit4->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,4))/(int)pow(2,4));
	       hMcWrongTrkBit5->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,5))/(int)pow(2,5));
	       hMcWrongTrkBit8->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,8))/(int)pow(2,8));  
	       hMcWrongTrkBit9->Fill((Reco_mu_trk_PIDmask[imutr] & (int)pow(2,9))/(int)pow(2,9));
	       } */ 
	  }
	}
        if (Reco_QQ_type[iqq] == 2) {
	  QQMass2Trk_passmu3->Fill(theMass);
          QQPt2Trk_passmu3->Fill(theQQ4mom->Perp());
          if ( Reco_QQ_mulpt[iqq] == theMCMatchedTrkMu || Reco_QQ_muhpt[iqq] == theMCMatchedTrkMu ) {
	    hMcRightTrkTrkMuMass->Fill(theMass);  
	  } else {
	    hMcWrongTrkTrkMuMass->Fill(theMass);
	  }
          if (iqq == bestQQ) {
	    QQMass2Trk_best->Fill(theMass);
	  }
	}
	if (Reco_QQ_type[iqq] == 3) {
	  QQMass1Glob1Cal_passmu3->Fill(theMass);
          QQPt1Glob1Cal_passmu3->Fill(theQQ4mom->Perp());
	  if (iqq == bestQQ) QQMass1Glob1Cal_best->Fill(theMass);
	  //
	  // MC matched calo mu (after cuts)
	  //
	  /* if (Reco_QQ_mulpt[iqq] == theMCMatchedCalMu && (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2) ) {
	  //
	  hMcRightCalGlobMuMass->Fill(theMass);
	  hMcRightCalGlobMuDeltaR->Fill(Reco_QQ_DeltaR[iqq]);
	  hMcRightCalGlobMuVtxChi2->Fill(Reco_QQ_normChi2[iqq]);
	  hMcRightCalGlobMuS->Fill(Reco_QQ_s[iqq]);
	  hMcRightCalGlobMucosAlpha->Fill(Reco_QQ_cosAlpha[iqq]);
	  } else {
	  hMcWrongCalGlobMuMass->Fill(theMass);
	  hMcWrongCalGlobMuDeltaR->Fill(Reco_QQ_DeltaR[iqq]);
	  hMcWrongCalGlobMuVtxChi2->Fill(Reco_QQ_normChi2[iqq]);
	  hMcWrongCalGlobMuS->Fill(Reco_QQ_s[iqq]);
	  hMcWrongCalGlobMucosAlpha->Fill(Reco_QQ_cosAlpha[iqq]);
	  } */
	}
	if (Reco_QQ_type[iqq] == 0 && theQQ4mom->Perp() < 6.0) QQMass2GlobPT6_passmu3->Fill(theMass);
	if (Reco_QQ_type[iqq] == 1 && theQQ4mom->Perp() < 6.0) QQMass1Glob1TrkPT6_passmu3->Fill(theMass);
	if (Reco_QQ_type[iqq] == 3 && theQQ4mom->Perp() < 6.0) QQMass1Glob1CalPT6_passmu3->Fill(theMass);
          
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
	//	if (HLTBits_accept[0]) {
	  if (Reco_QQ_type[iqq] == 0) WSMass2Glob_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 1) WSMass1Glob1Trk_passmu3->Fill(theMass);
	  if (Reco_QQ_type[iqq] == 3) WSMass1Glob1Cal_passmu3->Fill(theMass);
	  // }
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

  // trigger passed 
  QQMass2Glob_passmu3              = new TH1F("QQMass2Glob_passmu3",  "Invariant mass (2 global muons)", 20, 2. ,4.);
  QQMass2Trk_passmu3               = new TH1F("QQMass2Trk_passmu3",  "Invariant mass (2 tracker muons)", 20, 2. ,4.);
  QQMass1Glob1Trk_passmu3          = new TH1F("QQMass1Glob1Trk_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 20, 2., 4.);
  QQMass1Glob1Cal_passmu3          = new TH1F("QQMass1Glob1Cal_passmu3",  "Invariant mass (1 global + 1 calo muon)", 20, 2. ,4.);
  QQPt2Glob_passmu3                = new TH1F("QQPt2Glob_passmu3",  "Pt (2 global muons)", 20, 0.,20.);
  QQPt2Trk_passmu3                 = new TH1F("QQPt2Trk_passmu3",  "Pt (2 tracker muons)", 20, 0.,20.);
  QQPt1Glob1Trk_passmu3            = new TH1F("QQPt1Glob1Trk_passmu3",  "Pt (1 global + 1 tracker muon)", 20, 0.,20.);
  QQPt1Glob1Cal_passmu3            = new TH1F("QQPt1Glob1Cal_passmu3",  "Pt (1 global + 1 calo muon)", 20, 0.,20.); 
  QQMass2Glob_best              = new TH1F("QQMass2Glob_best",  "Invariant mass (2 global muons)", 20, 2.,4.);
  QQMass2Trk_best              = new TH1F("QQMass2Trk_best",  "Invariant mass (2 tracker muons)", 20, 2.,4.);
  QQMass1Glob1Trk_best          = new TH1F("QQMass1Glob1Trk_best",  "Invariant mass (1 global + 1 tracker muon)", 20, 2.,4.);
  QQMass1Glob1Cal_best          = new TH1F("QQMass1Glob1Cal_best",  "Invariant mass (1 global + 1 calo muon)", 20, 2.,4.);
  QQMass2GlobPT6_passmu3              = new TH1F("QQMass2GlobPT6_passmu3",  "Invariant mass (2 global muons)", 20, 2.,4.);
  QQMass1Glob1TrkPT6_passmu3          = new TH1F("QQMass1Glob1TrkPT6_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 20, 2.,4.);
  QQMass1Glob1CalPT6_passmu3          = new TH1F("QQMass1Glob1CalPT6_passmu3",  "Invariant mass (1 global + 1 calo muon)", 20, 2.,4.);
  WSMass2Glob_passmu3              = new TH1F("WSMass2Glob_passmu3",  "Invariant mass (2 global muons)", 20, 2.,4.);
  WSMass1Glob1Trk_passmu3          = new TH1F("WSMass1Glob1Trk_passmu3",  "Invariant mass (1 global + 1 tracker muon)", 20, 2.,4.);
  WSMass1Glob1Cal_passmu3          = new TH1F("WSMass1Glob1Cal_passmu3",  "Invariant mass (1 global + 1 calo muon)", 20, 2.,4.);
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

  // deltaR
  hMcRecoGlobMuDeltaR                  = new TH1F("hMcRecoGlobMuDeltaR",  "MC-reco matching #Delta R (global muons)", 100, 0.,0.5);
  hMcRecoTrkMuDeltaR                  = new TH1F("hMcRecoTrkMuDeltaR",  "MC-reco matching #Delta R (tracker muons)", 100, 0.,0.5);
  hMcRecoCalMuDeltaR                  = new TH1F("hMcRecoCalMuDeltaR",  "MC-reco matching #Delta R (calo muons)", 100, 0.,0.5);

  // mc-truth matching
  hMcRightGlbMunPixHits            = new TH1F("hMcRightGlbMunPixHits",  "number of pixel hits - MC matched (global muons)", 8, -0.5, 7.5);
  hMcWrongGlbMunPixHits            = new TH1F("hMcWrongGlbMunPixHits",  "number of pixel hits - MC unmatched (global muons)", 8, -0.5, 7.5);
  hMcRightGlbMud0            = new TH1F("hMcRightGlbMud0",  "d0 - MC matched (global muons)", 100, 0., 30.);
  hMcWrongGlbMud0            = new TH1F("hMcWrongGlbMud0",  "d0 - MC unmatched (global muons)", 100, 0., 30.);
  hMcRightGlbMudz            = new TH1F("hMcRightGlbMudz",  "dz - MC matched (global muons)", 100, 0., 50.);
  hMcWrongGlbMudz            = new TH1F("hMcWrongGlbMudz",  "dz - MC unmatched (global muons)", 100, 0., 50.);
  hMcRightGlbMuFirstLayer             = new TH1F("hMcRightGlbMuFirstLayer",  "first pixel layer hit - MC matched (global muons)", 4, -0.5, 3.5);
  hMcWrongGlbMuFirstLayer             = new TH1F("hMcWrongGlbMuFirstLayer",  "first pixel layer hit - MC unmatched (global muons)", 4, -0.5, 3.5);
  hMcRightGlbGlbMuVtxProb                  = new TH1F("hMcRightGlbGlbMuVtxProb",  "Vertex probability - MC matched (global+global muons)", 1000, 0.0,1.0);
  hMcWrongGlbGlbMuVtxProb                  = new TH1F("hMcWrongGlbGlbMuVtxProb",  "Vertex probability - MC unmatched (global+global muons)", 1000, 0.0,1.0);
  hMcRightGlbGlbMuMass                  = new TH1F("hMcRightGlbGlbMuMass",  "Inv. mass - MC matched (global+global muons)", 20, 2.0,4.0);
  hMcWrongGlbGlbMuMass                  = new TH1F("hMcWrongGlbGlbMuMass",  "Inv. mass - MC unmatched (global+global muons)", 20, 2.0,4.0);
  hMcRightTrkTrkMuMass                  = new TH1F("hMcRightTrkTrkMuMass",  "Inv. mass - MC matched (tracker+tracker muons)", 20, 2.0,4.0);
  hMcWrongTrkTrkMuMass                  = new TH1F("hMcWrongTrkTrkMuMass",  "Inv. mass - MC unmatched (tracker+tracker muons)", 20, 2.0,4.0);
 
  // mc truth matching - trk
  hMcRightTrkMuPt                      = new TH1F("hMcRightTrkMuPt",  "pT  - MC matched (tracker muons)", 50, 0.,5.0);
  hMcWrongTrkMuPt                      = new TH1F("hMcWrongTrkMuPt",  "pT - MC unmatched (tracker muons)", 50, 0.,5.0);
  hMcRightTrkMuChi2                    = new TH1F("hMcRightTrkMuChi2",  "chi2  - MC matched (tracker muons)", 50, -0.5, 6.5);
  hMcWrongTrkMuChi2                      = new TH1F("hMcWrongTrkMuChi2",  "chi2 - MC unmatched (tracker muons)", 50, -0.5,6.5);
  hMcRightTrkMuNhits                    = new TH1F("hMcRightTrkMuNhits",  "chi2  - MC matched (tracker muons)", 30, 0.5, 30.5);
  hMcWrongTrkMuNhits                    = new TH1F("hMcWrongTrkMuNhits",  "chi2 - MC unmatched (tracker muons)", 30, 0.5,30.5);
  hMcRightTrkBit4              = new TH1F("hMcRightTrkBit4",  "2DCompatibilityLoose bit - MC matched (tracker muons)", 4, -1.5,2.5);
  hMcWrongTrkBit4              = new TH1F("hMcWrongTrkBit4",  "2DCompatibilityLoose bit - MC unmatched (tracker muons)", 4, -1.5,2.5);
  hMcRightTrkBit5              = new TH1F("hMcRightTrkBit5",  "2DCompatibilityTight bit - MC matched (tracker muons)", 4, -1.5,2.5);
  hMcWrongTrkBit5              = new TH1F("hMcWrongTrkBit5",  "2DCompatibilityTight bit - MC unmatched (tracker muons)", 4, -1.5,2.5);
  hMcRightTrkBit8              = new TH1F("hMcRightTrkBit8",  "StationOptimizedLowPtLoose bit - MC matched (tracker muons)", 4, -1.5,2.5);
  hMcWrongTrkBit8              = new TH1F("hMcWrongTrkBit8",  "StationOptimizedLowPtLoose bit - MC unmatched (tracker muons)", 4, -1.5,2.5);
  hMcRightTrkBit9              = new TH1F("hMcRightTrkBit9",  "StationOptimizedLowPtTight bit - MC matched (tracker muons)", 4, -1.5,2.5);
  hMcWrongTrkBit9              = new TH1F("hMcWrongTrkBit9",  "StationOptimizedLowPtTight bit - MC unmatched (tracker muons)", 4, -1.5,2.5);
  hMcRightTrkMud0            = new TH1F("hMcRightTrkMud0",  "d0 - MC matched (global muons)", 100, 0., 30.);
  hMcWrongTrkMud0            = new TH1F("hMcWrongTrkMud0",  "d0 - MC unmatched (global muons)", 100, 0., 30.);
  hMcRightTrkMudz            = new TH1F("hMcRightTrkMudz",  "dz - MC matched (global muons)", 100, 0., 50.);
  hMcWrongTrkMudz            = new TH1F("hMcWrongTrkMudz",  "dz - MC unmatched (global muons)", 100, 0., 50.);
  hMcRightGlbTrkMuVtxProb                  = new TH1F("hMcRightGlbTrkMuVtxProb",  "Vertex probability - MC matched (global+tracker muons)", 1000, 0.0,1.0);
  hMcWrongGlbTrkMuVtxProb                   = new TH1F("hMcWrongGlbTrkMuVtxProb",  "Vertex probability - MC unmatched (global+tracker muons)", 1000, 0.0,1.0);
  hMcRightGlbTrkMuMass                  = new TH1F("hMcRightGlbTrkMuMass",  "Inv. mass - MC matched (global+global muons)", 20, 2.0,4.0);
  hMcWrongGlbTrkMuMass                  = new TH1F("hMcWrongGlbTrkMuMass",  "Inv. mass - MC unmatched (global+global muons)", 20, 2.0,4.0);

  // mc truth matching - calo
  hMcRightCalMuPt                      = new TH1F("hMcRightCalMuPt",  "pT  - MC matched (calo muons)", 50, 0.,5.0);
  hMcWrongCalMuPt                      = new TH1F("hMcWrongCalMuPt",  "pT - MC unmatched (calo muons)", 50, 0.,5.0);
  hMcRightCalMuCaloComp                      = new TH1F("hMcRightCalMuCaloComp",  "pT  - MC matched (calo muons)", 60, 0.5,1.1);
  hMcWrongCalMuCaloComp                      = new TH1F("hMcWrongCalMuCaloComp",  "pT - MC unmatched (calo muons)", 60, 0.5,1.1);
  hMcRightCalMuChi2                    = new TH1F("hMcRightCalMuChi2",  "chi2  - MC matched (calo muons)", 50, -0.5, 6.5);
  hMcWrongCalMuChi2                      = new TH1F("hMcWrongCalMuChi2",  "chi2 - MC unmatched (calo muons)", 50, -0.5,6.5);
  hMcRightCalMuNhits                    = new TH1F("hMcRightCalMuNhits",  "chi2  - MC matched (calo muons)", 30, 0.5, 30.5);
  hMcWrongCalMuNhits                      = new TH1F("hMcWrongCalMuNhits",  "chi2 - MC unmatched (calo muons)", 30, 0.5,30.5);
  hMcRightCalGlobMuDeltaR                  = new TH1F("hMcRightCalGlobMuDeltaR",  " DeltaR - MC matched (calo+global muons)", 80, 0.,5.0);
  hMcWrongCalGlobMuDeltaR                  = new TH1F("hMcWrongCalGlobMuDeltaR",  " DeltaR - MC unmatched (calo+global muons)", 80, 0.,5.0);
  hMcRightCalGlobMuMass                  = new TH1F("hMcRightCalGlobMuMass",  "Inv. mass - MC matched (calo+global muons)", 20, 2.0,4.0);
  hMcWrongCalGlobMuMass                  = new TH1F("hMcWrongCalGlobMuMass",  "Inv. Mass - MC unmatched (calo+global muons)", 20, 2.0,4.0);
  hMcRightCalGlobMuVtxChi2                  = new TH1F("hMcRightCalGlobMuVtxChi2",  "Vertex norm. #chi^{2} - MC matched (calo+global muons)", 80, 0.0,10.0);
  hMcWrongCalGlobMuVtxChi2                  = new TH1F("hMcWrongCalGlobMuVtxChi2",  "Vertex norm. #chi^{2} - MC unmatched (calo+global muons)", 80, 0.0,10.0);
  hMcRightCalGlobMuS                  = new TH1F("hMcRightCalGlobMuS",  "Significance of muon IPs - MC matched (calo+global muons)", 80, 0.0,100.0);
  hMcWrongCalGlobMuS                  = new TH1F("hMcWrongCalGlobMuS",  "Significance of muon IPs - MC unmatched (calo+global muons)", 80, 0.0,100.0);
  hMcRightCalGlobMucosAlpha                  = new TH1F("hMcRightCalGlobMucosAlpha",  "cos #alpha - MC matched (calo+global muons)", 100, -1.0,1.0);
  hMcWrongCalGlobMucosAlpha                  = new TH1F("hMcWrongCalGlobMucosAlpha",  "cos #alpha - MC unmatched (calo+global muons)", 100, -1.0,1.0);

  // all
  hMcRightAllMuIso                  = new TH1F("hMcRightAllMuIso",  "Isolation - MC matched (all muons)", 50, -1.0,15.0);
  
  // sizes 
  passingHLTMu3      = new TH1F("passingHLTMu3", "Pass HLTMu3", 4, -1.5,2.5);
  trackSize          = new TH1F("trackSize", "Track size", 100, 0.,100.);
  gammaSize          = new TH1F("gammaSize", "Gamma size", 100, 0.,100.);
  genQQSize          = new TH1F("genQQSize", "Generated J/psi size", 4, -1.5,2.5);
  genQQPt            = new TH1F("genQQPt", "Generated J/psi Pt", 10, 0.,20.);
  globalMuonSize     = new TH1F("globalmuonSize", "Global muon size", 5, -1.5,3.5);
  globalMuonPt       = new TH1F("globalmuonPt", "Global muon Pt", 10, 0.,20.);
  trackerMuonSize    = new TH1F("trackermuonSize", "Tracker muon size", 5, -1.5,3.5);
  trackerMuonPt      = new TH1F("trackermuonPt", "Tracker muon Pt", 10, 0.,10.);
  caloMuonSize       = new TH1F("calomuonSize", "Calo muon size", 5, -1.5,3.5); 
  caloMuonPt         = new TH1F("calomuonPt", "Calo muon Pt", 10, 0.,10.);
}

void finalJPsiAnalysis::drawPlots() {

  // TCanvas c; c.cd();
  // EleHistoGt4_size              -> Draw();  // c.Print("EleGt4.eps");
  
}

void finalJPsiAnalysis::saveHistos() {

  TFile fOut("Outfile_histo.root", "RECREATE");
  fOut.cd();
  
  QQMass2Glob_passmu3        -> Write();
  QQMass2Trk_passmu3         -> Write();
  QQMass1Glob1Trk_passmu3    -> Write(); 
  QQMass1Glob1Cal_passmu3    -> Write(); 
  QQPt2Glob_passmu3          -> Write(); 
  QQPt2Trk_passmu3           -> Write();
  QQPt1Glob1Trk_passmu3      -> Write(); 
  QQPt1Glob1Cal_passmu3      -> Write(); 
  QQMass2Glob_best           -> Write();
  QQMass2Trk_best           -> Write();
  QQMass1Glob1Trk_best       -> Write(); 
  QQMass1Glob1Cal_best       -> Write();
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
  hMcRecoGlobMuDeltaR       -> Write();   
  hMcRecoTrkMuDeltaR        -> Write();   
  hMcRecoCalMuDeltaR        -> Write();
  hMcRightGlbMunPixHits   -> Write(); 
  hMcWrongGlbMunPixHits   -> Write(); 
  hMcRightGlbMud0         -> Write(); 
  hMcWrongGlbMud0         -> Write(); 
  hMcRightGlbMudz         -> Write(); 
  hMcWrongGlbMudz         -> Write(); 
  hMcRightGlbMuFirstLayer -> Write(); 
  hMcWrongGlbMuFirstLayer -> Write();
  hMcRightGlbTrkMuVtxProb -> Write();   
  hMcWrongGlbTrkMuVtxProb -> Write(); 
  hMcRightGlbTrkMuMass -> Write();   
  hMcWrongGlbTrkMuMass -> Write(); 
  hMcRightGlbGlbMuVtxProb -> Write(); 
  hMcWrongGlbGlbMuVtxProb -> Write();
  hMcRightGlbGlbMuMass -> Write(); 
  hMcWrongGlbGlbMuMass -> Write();
  hMcRightTrkTrkMuMass -> Write(); 
  hMcWrongTrkTrkMuMass -> Write();
  hMcRightTrkMuPt   -> Write(); 
  hMcWrongTrkMuPt   -> Write(); 
  hMcRightTrkBit4   -> Write(); 
  hMcWrongTrkBit4   -> Write(); 
  hMcRightTrkBit5   -> Write(); 
  hMcWrongTrkBit5   -> Write(); 
  hMcRightTrkBit8   -> Write(); 
  hMcWrongTrkBit8   -> Write(); 
  hMcRightTrkBit9   -> Write(); 
  hMcWrongTrkBit9    -> Write();
  hMcRightTrkMuChi2  -> Write();   
  hMcWrongTrkMuChi2  -> Write(); 
  hMcRightTrkMuNhits -> Write(); 
  hMcWrongTrkMuNhits -> Write(); 
  hMcRightCalMuChi2  -> Write(); 
  hMcWrongCalMuChi2  -> Write(); 
  hMcRightCalMuNhits -> Write(); 
  hMcWrongCalMuNhits -> Write(); 
  hMcRightTrkMud0         -> Write(); 
  hMcWrongTrkMud0         -> Write(); 
  hMcRightTrkMudz         -> Write(); 
  hMcWrongTrkMudz         -> Write(); 
  hMcRightCalMuCaloComp     -> Write();   
  hMcWrongCalMuCaloComp     -> Write(); 
  hMcRightCalMuPt           -> Write(); 
  hMcWrongCalMuPt           -> Write();
  hMcRightAllMuIso          -> Write();
  hMcRightCalGlobMuDeltaR   -> Write();   
  hMcWrongCalGlobMuDeltaR   -> Write();   
  hMcRightCalGlobMuMass     -> Write();   
  hMcWrongCalGlobMuMass     -> Write(); 
  hMcRightCalGlobMuVtxChi2  -> Write(); 
  hMcWrongCalGlobMuVtxChi2  -> Write();
  hMcRightCalGlobMuS  -> Write(); 
  hMcWrongCalGlobMuS  -> Write();
  hMcRightCalGlobMucosAlpha  -> Write(); 
  hMcWrongCalGlobMucosAlpha  -> Write();

  passingHLTMu3     -> Write();
  trackSize         -> Write();
  gammaSize         -> Write();
  genQQSize         -> Write();
  genQQPt           -> Write();
  globalMuonSize    -> Write();
  globalMuonPt      -> Write();
  trackerMuonSize   -> Write();
  trackerMuonPt     -> Write();
  caloMuonSize      -> Write();
  caloMuonPt        -> Write();
}			      
			       
double finalJPsiAnalysis::PhiInRange(const double& phi) const {
      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}

double finalJPsiAnalysis::deltaR(const TLorentzVector* t, const TLorentzVector* u) const {
      return sqrt(pow(t->Eta()-u->Eta(),2) +pow(PhiInRange(t->Phi()-u->Phi()),2));
}

int finalJPsiAnalysis::theBestQQ() {
    
 int theBest = -1;
  float thehighestPt = -1.;
 
  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 0 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];   if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];   if (thelptMu >= Reco_mu_glb_size) continue;
      if (Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi && 
	  Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb && 
	  (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) && 
	  fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	  Reco_mu_glb_nhitstrack[thelptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thelptMu] < MAX_normchi2_glb &&
	  (((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thelptMu] == 1)) &&
	  fabs(Reco_mu_glb_d0[thelptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thelptMu]) < MAX_dz_trk
	  ) {
	return iqq;
      }
    }
  }

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 1 ) {
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;

      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
	   Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) &&
	   fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
        if (theTrMumom->Perp() > thehighestPt) {
	  thehighestPt = theTrMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 2 ) {
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_trk_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;

      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
	   Reco_mu_trk_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thehptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thehptMu] + Reco_mu_trk_nhitsPixE[thehptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thehptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thehptMu]) < MAX_dz_trk &&
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thehptMu);
        if (theTrMumom->Perp() > thehighestPt) {
	  thehighestPt = theTrMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 3 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];
      int thelptMu = Reco_QQ_mulpt[iqq];
      if (thelptMu >= Reco_mu_cal_size) {
	// cout << "Non deve succedere! cmIndex = " << thelptMu+1 << " cmSize = " << Reco_mu_cal_size << endl;
	continue;
      }

      if ( Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   Reco_mu_cal_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   Reco_mu_cal_normChi2[thelptMu] < 3.0 && 
	   Reco_mu_cal_caloComp[thelptMu] > 0.89) {
	
        TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(thelptMu);
        if (theCaMumom->Perp() > thehighestPt) {
	  thehighestPt = theCaMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  return theBest;
} 
	
	 
