/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
  * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * > root -q -b 
 * > .L RooModifTsallis.cc+ 
 * > .x generateComplicatedPtTemplates.C+
 * 2D templates are written to "destDir"
 *
 */

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include <sstream>
#include <vector>
#include <iostream>
#include <stdio.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
//#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif

#include "RooRealVar.h"
// #include "RooModifTsallis.h"
//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//---
int useSqrts=0;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TString melaName = "ZZLD"; // name of MELA branch to be used. Possibilities are ZZLD,ZZLD_analBkg,ZZLD_postICHEP,ZZLD_PtY,pseudoMelaLD, spinTwoMinimalMelaLD 

bool extendToHighMass = true; // Include signal samples above 600 GeV

float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;

const TString destDir = "../CreateDatacards/templates2D/";
// static const int nsamp = 8;
// const TString dataFileNames[nsamp] = {"gg","vbf","wh","zh","tth","zz","zx","ggzz"};
static const int nsamp = 5;
const TString dataFileNames[nsamp] = {"gg","vbf","wh","zh","zz"};
TString systSources[nsamp][5];

//=======================================================================

TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,int LHCsqrts= 7, TString systName = "Default", bool down = false){  
 
  char fileName[200];

  TFile* f8compl = new TFile(destDir + "PtoverMComplicated_mZZ_" + dataFileNames[sampleIndex] + "_" + channel + "_8TeV.root");

  TFile* f7simpl = new TFile(destDir + "PtoverM_mZZ_" + dataFileNames[sampleIndex] + "_" + channel + "_7TeV.root");

  if (systName == "OneSyst") {
    if (down) sprintf(fileName,"h_Ptmzz_mzz_OneSyst_down");
    else sprintf(fileName,"h_Ptmzz_mzz_OneSyst_up");
  } else sprintf(fileName,"h_Ptmzz_mzz");
  cout << fileName << endl;
  
  TH2F* Hist8c = (TH2F*)f8compl->Get(fileName);
  TH2F* Hist7s = (TH2F*)f7simpl->Get(fileName);

  TH2F* bkgHist = (TH2F*)Hist8c->Clone();
  
  double norm;
  TH1F* tempProj;
  TH1F* tempProj7;
  tempProj = (TH1F*)Hist8c->ProjectionY("tempProj",13,13);   // 125-126 GeV
  tempProj7 = (TH1F*)Hist7s->ProjectionY("tempProj7",13,13);   // 125-126 GeV
  TH1F* weight = (TH1F*)tempProj7->Clone();
  weight->Divide(tempProj7,tempProj);

  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();

  for(int i=1; i<=nXbins; i++){
    for(int j=1; j<=nYbins; j++){
      bkgHist->SetBinContent(i,j,bkgHist->GetBinContent(i,j)*weight->GetBinContent(j));
    }
  }

  // normalize slices
 
  for(int i=1; i<=nXbins; i++){
    
    tempProj = (TH1F*)bkgHist->ProjectionY("tempProj",i,i);
    norm = tempProj->Integral();
    // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << norm << endl;

    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	bkgHist->SetBinContent(i,j,bkgHist->GetBinContent(i,j)/norm);
      }
    }

  }
  
  //f8compl->Close();
  //f7simpl->Close();

  return bkgHist;
  
}

//=======================================================================

void makeTemplate(TString channel="4mu"){

  char fileName[200];

  for (Int_t lhc=7; lhc<8; lhc++) {
    for (Int_t k=0; k<nsamp; k++) {
     
      TString lhcs = "7";
      if (lhc==8) lhcs="8";

      TFile* ftemp = new TFile(destDir + "PtoverMComplicated_mZZ_" + dataFileNames[k] + "_" + channel + "_"+ lhcs + "TeV.root","RECREATE");

      cout << "*** Now filling: " << dataFileNames[k] << ", Default template" << endl;
      TH2F* h_mzzD = fillTemplate(channel,k,lhc);
     
      ftemp->cd();
  
      h_mzzD->Write("h_Ptmzz_mzz");
      cout << " " << h_mzzD->GetBinContent(100,10) << endl;
      
      cout << "*** Now filling: " << dataFileNames[k] << ", OneSyst templates" << endl;
      
      TH2F* h_mzzDup = fillTemplate(channel,k,lhc,"OneSyst",false);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_up");
      ftemp->cd();
      h_mzzDup->Write(fileName);
      cout << " " << h_mzzDup->GetBinContent(100,10) << endl;
      
      TH2F* h_mzzDdown = fillTemplate(channel,k,lhc,"OneSyst",true);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_down");
      ftemp->cd();
      h_mzzDdown->Write(fileName);
      cout << " " << h_mzzDdown->GetBinContent(100,10) << endl;

      ftemp->Close();

    }
  }

}

//=======================================================================

void generateComplicatedPtTemplates7TeV(){

  // bool isLowMass = false;

  /* systSources[0][0] = "Resummation";
  systSources[0][1] = "TopMass";
  systSources[0][2] = "Mela";
  
  systSources[1][0] = "PDF-VBF";
  systSources[1][1] = "scale-VBF";
  systSources[1][2] = "Mela";
  
  systSources[5][0] = "SingleZ";
  systSources[5][1] = "PDF-ZZ";
  systSources[5][2] = "scale-ZZ";
  systSources[5][3] = "Mela"; */

  makeTemplate("4mu");
  makeTemplate("4e");
  makeTemplate("2e2mu");

}



