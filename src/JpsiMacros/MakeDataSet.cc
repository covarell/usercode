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
  const float JpsiPtMin = 1.;
  const float JpsiPtMax = 100.;

  // RooFit stuff
  RooRealVar* JpsiMass = new RooRealVar("JpsiMass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar* MCweight = new RooRealVar("MCweight","Monte Carlo Weight",0.,5.);

  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi Pt",JpsiPtMin,JpsiPtMax,"GeV/c");

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

  RooDataSet* data = new RooDataSet("data","Prompt sample",RooArgList(*JpsiMass,*Jpsict,*MCweight,JpsiType,MCType,*JpsiPt));

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
    else if(MCcat == 1) weight = 0.1745;
    else if(MCcat == 2) weight = 2.2831;

    //exclude real Jpsi in background MC
    if(MCcat == 2){
      bool aRealJpsiEvent = false;
      for (int imugen=0; imugen<Mc_mu_size; imugen++) {
	if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
      }
      if (aRealJpsiEvent) continue;
    }

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

      if (Reco_QQ_sign[iqq] != 0) continue;

      if(Reco_QQ_type[iqq] == 2 || Reco_QQ_type[iqq] > 3) continue;

	TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
	float theMass = theQQ4mom->M();
        float theCtau = Reco_QQ_ctau[iqq]*10.;
	float thePt = theQQ4mom->Pt();

        if (theMass > JpsiMassMin && theMass < JpsiMassMax && theCtau > JpsiCtMin && theCtau < JpsiCtMax && thePt > JpsiPtMin && thePt < JpsiPtMax) {

	  passedCandidates++;

	  JpsiMass->setVal(theMass);
	  Jpsict->setVal(Reco_QQ_ctau[iqq]*10.);
	  JpsiType.setIndex(Reco_QQ_type[iqq],kTRUE);
	  JpsiPt->setVal(thePt);

	  MCType.setIndex(MCcat,kTRUE);
	  MCweight->setVal(weight);

	  data->add(RooArgSet(*JpsiMass,*Jpsict,*MCweight,JpsiType,MCType,*JpsiPt));

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

  c1.SaveAs("firstTest.eps");

  TFile fOut("Out_DataSet.root", "RECREATE");
  fOut.cd();
  data->Write();
  fOut.Close();

  // summary
  cout << "total number of events = " << totalEvents << endl;
  cout << "number of passed candidates = " << passedCandidates << endl;

} // end of program
			       
			       
