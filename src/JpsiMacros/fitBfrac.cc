#include <TStyle.h>
#include <TCanvas.h>
#include "TLine.h"
#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"

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

#include "fitBfrac.hh"

using namespace std;
using namespace RooFit;

fitBfrac::fitBfrac(TTree *tree) 
  : JPsiTreeBase(tree) {
}

fitBfrac::~fitBfrac(){ } 

void fitBfrac::Loop() {

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries();
  
  // booking histos
  bookHistos();

  // loop over events
  std::cout << "Number of entries = " << nentries << std::endl;

  // counters
  int totalEvents = 0;

  // mass-lifetime limits
  float JpsiMassMin = 2.6;
  float JpsiMassMax = 3.5;
  float JpsiCtMin = -1.0;
  float JpsiCtMax = 5.0;

  // RooFit stuff
  RooRealVar* JpsiMass = new RooRealVar("JpsiMass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi pt",0.,60.,"GeV/c");
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooDataSet* dataGGNonP = new RooDataSet("dataGGNonP",
					  "Non-prompt sample (glob-glob)",
					  RooArgList(*JpsiMass,*Jpsict));
  RooDataSet* dataGTNonP = new RooDataSet("dataGTNonP",
					  "Non-prompt sample (glob-trk)",
					  RooArgList(*JpsiMass,*Jpsict));
  RooDataSet* dataGCNonP = new RooDataSet("dataGCNonP",
					  "Non-prompt sample (glob-calo)",
					  RooArgList(*JpsiMass,*Jpsict));
                
  for (int jentry=0; jentry<nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    totalEvents++;

    // bool aRealJpsiEvent = false;
    // for (int imugen=0; imugen<Mc_mu_size; imugen++) {
    //   if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
    // }
    // if (aRealJpsiEvent) continue;

    // TRIGGER CUTS 
    if (!HLTBits_accept[0]) continue;    // SingleMu3

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {
      if (Reco_QQ_sign[iqq] == 0) {
	TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
	float theMass = theQQ4mom->M();
        float theCtau = Reco_QQ_ctau[iqq]*10.;
        if (theMass > JpsiMassMin && theMass < JpsiMassMax && theCtau > JpsiCtMin && theCtau < JpsiCtMax) {
	  if (Reco_QQ_type[iqq] == 0) {
	    JpsiMass->setVal(theMass);
	    Jpsict->setVal(Reco_QQ_ctau[iqq]*10.);
	    dataGGNonP->add(RooArgSet(*JpsiMass,*Jpsict));
	  }
          if (Reco_QQ_type[iqq] == 1) {
	    JpsiMass->setVal(theMass);
	    Jpsict->setVal(Reco_QQ_ctau[iqq]*10.);
	    dataGTNonP->add(RooArgSet(*JpsiMass,*Jpsict));
	  }
          if (Reco_QQ_type[iqq] == 3) {
	    JpsiMass->setVal(theMass);
	    Jpsict->setVal(Reco_QQ_ctau[iqq]*10.);
	    dataGCNonP->add(RooArgSet(*JpsiMass,*Jpsict));
	  }
	}
      }
    }
  }

  TCanvas c1("c1","c1",10,10,600,1000);
  c1.Divide(2,3);

  c1.cd(1);
  RooPlot* frameMass = JpsiMass->frame();
  dataGGNonP->plotOn(frameMass,Binning(50));
  frameMass->Draw();

  c1.cd(3);
  RooPlot* frameMass2 = JpsiMass->frame();
  dataGTNonP->plotOn(frameMass2,Binning(50));
  frameMass2->Draw();

  c1.cd(5);
  RooPlot* frameMass3 = JpsiMass->frame(); 
  dataGCNonP->plotOn(frameMass3,Binning(50));
  frameMass3->Draw();

  c1.cd(2);
  gPad->SetLogy(1);
  RooPlot* framect = Jpsict->frame();
  dataGGNonP->plotOn(framect,Binning(50)); 
  framect->Draw();

  c1.cd(4);
  gPad->SetLogy(1);
  RooPlot* framect2 = Jpsict->frame();
  dataGTNonP->plotOn(framect2,Binning(50)); 
  framect2->Draw();

  c1.cd(6);
  gPad->SetLogy(1);
  RooPlot* framect3 = Jpsict->frame();
  dataGCNonP->plotOn(framect3,Binning(50)); 
  framect3->Draw(); 

  c1.SaveAs("firstTest.eps");

  saveHistos();

  drawPlots();

  // summary
  cout << "total number of events = "              << totalEvents     << endl;

} // end of program


void fitBfrac::bookHistos() {

  // QQMass1Glob1Cal_pass2mu3ups          = new TH1F("QQMass1Glob1Cal_pass2mu3ups",  "Invariant mass (1 global + 1 calo muon)", 100, 0.,15.);
}

void fitBfrac::drawPlots() {

  // TCanvas c; c.cd();
  // EleHistoGt4_size              -> Draw();  // c.Print("EleGt4.eps");
  
}

void fitBfrac::saveHistos() {

  // TFile fOut("Outfile_histo.root", "RECREATE");
  // fOut.cd();
    
  // QQMass1Glob1Cal_pass2mu3ups  -> Write();  
}			      
			       
			       
