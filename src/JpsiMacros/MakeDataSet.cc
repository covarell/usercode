#include <TStyle.h>
#include <TCanvas.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"

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

#include "MakeDataSet.hh"

using namespace std;
using namespace RooFit;

MakeDataSet::MakeDataSet(TTree *tree) 
  : JPsiTreeBase(tree) {
}

MakeDataSet::~MakeDataSet(){ } 

void MakeDataSet::Loop() {

  // to be defined as parameter in JPsiFitApp !!!
  bool onlyTheBest = true;

  MIN_nhits_trk = 12;
  MAX_normchi2_trk = 1.9;
  MAX_normchi2_glb = 10;

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // loop over events
  cout << "Number of entries = " << nentries << endl;

  // counters
  int totalEvents = 0;
  int passedCandidates = 0;

  // mass-lifetime limits
  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;
  const float JpsiCtMin = -1.0;
  const float JpsiCtMax = 5.0;

  // RooFit stuff
  RooRealVar* JpsiMass = new RooRealVar("JpsiMass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi pt",0.,60.,"GeV/c");
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar* MCweight = new RooRealVar("MCweight","Monte Carlo Weight",0.,5.);

  //categories for reconstruction
  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  //categories for MC production
  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  int MCcat = -999;

  float weight = 0.;

  RooDataSet* data = new RooDataSet("data","Prompt sample",RooArgList(*JpsiMass,*Jpsict,*MCweight,JpsiType,MCType));

  for (int jentry=0; jentry< nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    if (jentry%100000 == 0) cout << ">>> Processing event # " << jentry << endl;
    
    totalEvents++;

    // bool aRealJpsiEvent = false;
    // for (int imugen=0; imugen<Mc_mu_size; imugen++) {
    //   if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
    // }
    // if (aRealJpsiEvent) continue;

    // TRIGGER CUTS 
    if (!HLTBits_accept[0]) continue;    // SingleMu3

    //estrablish which kind of MC this is
    MCcat = -999;
    TString filestring(fChain->GetCurrentFile()->GetName());
    if(filestring.Contains("promptJpsiMuMu")) MCcat = 0;
    else if(filestring.Contains("inclBtoJpsiMuMu")) MCcat = 1;
    else if(filestring.Contains("InclusiveppToMu")) MCcat = 2;

    //set the MC weight for the different categories
    if(MCcat == 0) weight = 0.862;
    else if(MCcat == 1) weight = 0.1745*2;   // take into account b+bbar!!!
    else if(MCcat == 2) weight = 2.2831;
  
    //exclude real Jpsi in background MC
    if(MCcat == 2){
      bool aRealJpsiEvent = false;
      for (int imugen=0; imugen<Mc_mu_size; imugen++) {
	if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
      }
      if (aRealJpsiEvent) continue;
    }

    int theMCMatchedGlbMu1 = -1;
    int theMCMatchedGlbMu2 = -1;
    int theMCMatchedTrkMu = -1;
    int theMCMatchedCalMu = -1;
    unsigned int nGlobMuMatched = 0;

    //
    // Find MC matched muons (all QQ categories)
    //
    for (int imugen=0; imugen<Mc_mu_size; imugen++) {
      // cout << "Mc_mu_size = " << Mc_mu_size << endl;
      if (Mc_mumoth_id[imugen] == 443) {
        TLorentzVector *theMcMumom = (TLorentzVector*)Mc_mu_4mom->At(imugen);
        for (int imugl=0; imugl<Reco_mu_glb_size; imugl++) {
          // cout << "Reco_mu_glb_size = " << Reco_mu_glb_size << endl;
          TLorentzVector *theGlMumom = (TLorentzVector*)Reco_mu_glb_4mom->At(imugl);
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
          if (deltaR(theMcMumom,theTrMumom) < 0.03) {
	    theMCMatchedTrkMu = imutr;      break;
	  }
	}
         for (int imuca=0; imuca<Reco_mu_cal_size; imuca++) {
	   // cout << "Reco_mu_cal_size = " << Reco_mu_cal_size << endl;
          TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(imuca);
          if (deltaR(theMcMumom,theCaMumom) < 0.03) {
	    theMCMatchedCalMu = imuca;       break;
	  }
	}
      }      
    }

    // Find the best candidate (if needed)
    int myBest = 0;
    if (onlyTheBest) myBest = theBestQQ();

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

      if (Reco_QQ_sign[iqq] != 0) continue;

      if(Reco_QQ_type[iqq] == 2 || Reco_QQ_type[iqq] > 3) continue;

      if (onlyTheBest && iqq != myBest) continue;

	TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
	float theMass = theQQ4mom->M();
        float theCtau = Reco_QQ_ctau[iqq]*10.;

        if (theMass > JpsiMassMin && theMass < JpsiMassMax && theCtau > JpsiCtMin && theCtau < JpsiCtMax) {

	  passedCandidates++;

	  JpsiMass->setVal(theMass);
	  Jpsict->setVal(Reco_QQ_ctau[iqq]*10.);
	  JpsiType.setIndex(Reco_QQ_type[iqq],kTRUE);

          // Now, AFTER setting the weight, change to consider MC truth!
	  if (filestring.Contains("promptJpsiMuMu") || filestring.Contains("inclBtoJpsiMuMu")) {
	    bool isMatchedGlbGlb = (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu2) || (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2 && Reco_QQ_mulpt[iqq] == theMCMatchedGlbMu1);
            bool isMatchedGlbTrk = (Reco_QQ_mulpt[iqq] == theMCMatchedTrkMu && (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2) );
	    bool isMatchedGlbCal = (Reco_QQ_mulpt[iqq] == theMCMatchedCalMu && (Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu1 || Reco_QQ_muhpt[iqq] == theMCMatchedGlbMu2) ); 
	    if (!isMatchedGlbGlb && !isMatchedGlbTrk && !isMatchedGlbCal) MCcat = 2;
	  }

	  MCType.setIndex(MCcat,kTRUE);
	  MCweight->setVal(weight);

	  data->add(RooArgSet(*JpsiMass,*Jpsict,*MCweight,JpsiType,MCType));

	}
      }
  }

  data->setWeightVar(*MCweight);

  TCanvas c1("c1","c1",10,10,600,1000);
  c1.Divide(2,3);

  c1.cd(1);
  RooPlot* frameMass = JpsiMass->frame();
  data->plotOn(frameMass,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"));
  frameMass->Draw();

  c1.cd(3);
  RooPlot* frameMass2 = JpsiMass->frame();
  data->plotOn(frameMass2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"));
  frameMass2->Draw();

  c1.cd(5);
  RooPlot* frameMass3 = JpsiMass->frame(); 
  data->plotOn(frameMass3,Binning(50),RooFit::Cut("JpsiType==JpsiType::GC"));
  frameMass3->Draw();

  c1.cd(2);
  gPad->SetLogy(1);
  RooPlot* framect = Jpsict->frame();
  data->plotOn(framect,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"));
  framect->Draw();

  c1.cd(4);
  gPad->SetLogy(1);
  RooPlot* framect2 = Jpsict->frame();
  data->plotOn(framect2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"));
  framect2->Draw();

  c1.cd(6);
  gPad->SetLogy(1);
  RooPlot* framect3 = Jpsict->frame();
  data->plotOn(framect3,Binning(50),RooFit::Cut("JpsiType==JpsiType::GC"));
  framect3->Draw(); 

  c1.SaveAs("bestCands.gif");

  TFile fOut("DataSet_bestCands.root", "RECREATE");
  fOut.cd();
  data->Write();
  fOut.Close();

  // summary
  cout << "total number of events = " << totalEvents << endl;
  cout << "number of passed candidates = " << passedCandidates << endl;

} // end of program
		

	       
int MakeDataSet::theBestQQ() {
    
  int theBest = -1;
  float thehighestPt = -1.;
 
  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {
    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 0 &&
	Reco_mu_glb_nhitstrack[iqq] > MIN_nhits_trk && Reco_mu_glb_normChi2[iqq] < MAX_normchi2_glb) return iqq;
  }

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {
    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 1 ) {
      
      int theTM = Reco_QQ_mulpt[iqq];
      if (theTM >= Reco_mu_trk_size) {
	// cout << "Non deve succedere! tmIndex = " << theTM+1 << " tmSize = " << Reco_mu_trk_size << endl;
	continue;
      }

      if ( Reco_mu_trk_nhitstrack[theTM] > MIN_nhits_trk && ((Reco_mu_trk_PIDmask[theTM] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[theTM] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   Reco_mu_trk_normChi2[theTM] < MAX_normchi2_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(theTM);
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
      
      int theCM = Reco_QQ_mulpt[iqq];
      if (theCM >= Reco_mu_cal_size) {
	// cout << "Non deve succedere! cmIndex = " << theCM+1 << " cmSize = " << Reco_mu_cal_size << endl;
	continue;
      }

      if ( Reco_mu_cal_nhitstrack[theCM] > MIN_nhits_trk && Reco_mu_cal_normChi2[theCM] < 3.0 && Reco_mu_cal_caloComp[theCM] > 0.89) {
	
        TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(theCM);
        if (theCaMumom->Perp() > thehighestPt) {
	  thehighestPt = theCaMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  return theBest;
} 			       

double MakeDataSet::PhiInRange(const double& phi) const {
      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}

double MakeDataSet::deltaR(const TLorentzVector* t, const TLorentzVector* u) const {
      return sqrt(pow(t->Eta()-u->Eta(),2) +pow(PhiInRange(t->Phi()-u->Phi()),2));
}
