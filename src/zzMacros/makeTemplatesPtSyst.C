// Directly produce MC templates from samples: to be used with new MC

// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>

using namespace std;

static const int massRanges = 4;
static const int melaRanges = 3;

// static const int nWeightResum = 6;
// int massResum[nWeightResum] = {125,200,400,600,800,1000};
static const int nWeightResum = 2;
int massResum[nWeightResum] = {125,400};

static const int nWeightNLOvh = 3;
int massNLOvh[nWeightNLOvh] = {115,125,140};

float ptVar(float pt, float m, int overM = 0) {
  /* if (overM) return pt/m;
     return pt;	*/
  if (overM == 0) return pt;
  else if (overM == 1) return pt/m;
  else if (log(pt/m) < -5.5) return -5.4999;
  else if (log(pt/m) > 1.5) return 1.4999;
  return log(pt/m);
}

bool notVBFtagged(float njets) {
  if (njets >= 2) return false;
  return true;
}

void evalBinMigration(TH2F* def, TH2F* var, char fileName[200], bool syst = true, int binlimit = 20) { // i.e. 50 GeV
  
  char outpStr[200];
  sprintf(outpStr,"MCstatistics");
  // TH1F* def = (TH1F*)theFile->Get("ptH_Default");
  if (!def) return;
  float theErr;
  if (syst) {
    theErr = (def->Integral(binlimit,160)/def->Integral() - var->Integral(binlimit,160)/var->Integral())*(def->Integral()/def->Integral(binlimit,160));
    sprintf(outpStr,"%s",fileName);
  } else {
    theErr = 1./sqrt(def->GetEntries()*def->Integral(binlimit,160)/def->Integral());
  }
  cout << outpStr << " : " << int(theErr*10000)/100. << "%" << endl;
  return;
}

void adjustHistogram(TH2F* hist) {
  
  // transform in prob. dens.  
  /* for (Int_t iBin = 1; iBin <= hist->GetNbinsX() ; ++iBin) {
    for (Int_t iBiny = 1; iBiny <= hist->GetNbinsY() ; ++iBiny) {
      float a = hist->GetXaxis()->GetBinWidth(iBin)*hist->GetYaxis()->GetBinWidth(iBiny);
      hist->SetBinContent(iBin,iBiny,hist->GetBinContent(iBin,iBiny)/a);
      // if (hist->GetBinContent(iBin,iBiny) <= 0 ) hist->SetBinContent(iBin,iBiny,0.000001);
    }
    } */

  // normalize slices
  double norm;
  TH1F* tempProj;
  
  for(int i=1; i<=hist->GetNbinsX(); i++){
    
    tempProj = (TH1F*)hist->ProjectionY("tempProj",i,i);
    norm = tempProj->Integral();
 
    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=hist->GetNbinsX(); j++){
	hist->SetBinContent(i,j,hist->GetBinContent(i,j)/norm);
      }
    }

  }
  return;
}

int returnClosestMass(const int n, int* theArray, double mass) {

  double theMin = 999999.;
  int whichOne = -1;
  for (Int_t i = 1; i < n ; ++i) {
    if (fabs(mass - (float)theArray[i]) < theMin) {
      theMin = fabs(mass - (float)theArray[i]);
      whichOne = i;
    }
  }
  return whichOne;
}

void makeTemplatesPtSyst(char* = "2mu2e", int whichtype = 1, 
			int overM = 0, 
			bool also7TeV = true, bool moreFiles = true) { 
    
  cout << "overM = " << overM << endl;
     
  // SYSTEMATIC CODES            
  // -7 - NLO/LO : ZH
  // -6 - NLO/LO : WH
  // -5 - PDF : VBF
  // -4 - scales : VBF  
  // -3 - effect of finite top mass
  // -2 - mu_Q in gg signal shape
  // -1 - fraction of VBF (NOT for VBF fraction fitting)
  // 0 - DEFAULTS: ALWAYS RUN FIRST
  // 1 - Z+X in background
  // 2 - Cross-check Z+X vs. Z+jets
  // 3 - low pT shape in background: nonblinded
  // 4 - low pT shape in background: single Z
  // 5 - adding/not adding ggZZ
  // 6 - PDF : ZZ
  // 7 - scales : ZZ

  // VARIABLE BINNING (DO NOT USE)
  /* static const int ptmbins = 20;
   double ptmbinLimits[ptmbins+1] = {0.,.05,.1,.2,.3,.4,.5,.6,.7,
                                     .8,.9,1.0,1.1,1.2,1.3,1.4,1.5,
                                      1.75,2.0,2.25,2.5};

  if (!overM) {
    for (Int_t j = 0; j < ptmbins; j++) { 
      ptmbinLimits[j] *= 125.;
    }
  }
 
  static const int mbins = 102;
  double mbinLimits[mbins+1];

  for (Int_t i1 = 0; i1 < 41; i1++) { // 100-180 in 2 GeV steps
    mbinLimits[i1] = 100.+i1*2.;
  }
  for (Int_t i2 = 41; i2 < 83; i2++) {  // 180-600 in 10 GeV steps
    mbinLimits[i2] = 190.+(i2-41)*10.;
  } 
  for (Int_t i3 = 83; i3 < mbins+1; i3++) {  // 600-XXX in 50 GeV steps
    mbinLimits[i3] = 650.+(i3-83)*50.;
    } */ 

  // COSTANT BINNING (DEFAULT WITH LOG) 

  static const int ptmbins = 20;
  double ptmbinLimits[ptmbins+1] = {-5.5,-5.15,-4.8,-4.45,-4.1,
                                    -3.75,-3.4,
				    -3.05,-2.7,-2.35,-2.,
				    -1.65,-1.3,-0.95,-0.6,
				    -0.25,0.1,0.45,0.8,
				    1.15,1.5};
                                    // constant binning from -5.5 to 1.5
  
  if (overM == 0) {
    for (Int_t j = 0; j < ptmbins+1; j++) { 
      ptmbinLimits[j] = 41.0*(ptmbinLimits[j] + 5.5);
    }
  }
  else if (overM == 1) {
    for (Int_t j = 0; j < ptmbins+1; j++) { 
      ptmbinLimits[j] = 0.325*(ptmbinLimits[j] + 5.5);
    }
  }
 
  static const int mbins = 750;
  double mbinLimits[mbins+1];

  for (Int_t i1 = 0; i1 < 751; i1++) { // all in 2 GeV steps
    mbinLimits[i1] = 100.+i1*2.;
  }

  bool withNLOMela = false;
  char nameFile[200] = "SelectedTree";
  char nameFile2[200] = "";
  if (withNLOMela) sprintf(nameFile,"angles");

  int sqrts = 8;
  if (also7TeV) sqrts = 7;

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
 
  // READ TREES
  TChain* ggTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    if (moreFiles) {
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H210.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H210.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H210.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
    }
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H125.root");  ggTree->Add(nameFile2);
    if (moreFiles) {      
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H145.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H116.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H128.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H117.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H127.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H500.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H135.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H375.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H119.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H129.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H118.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H121.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H122.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H850.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H115.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H123.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H145.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H116.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H128.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H117.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H127.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H500.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H135.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H375.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H119.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H129.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H118.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H121.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H122.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H850.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H115.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H123.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H145.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H126.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H124.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H550.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H700.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H450.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H116.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H130.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H350.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H300.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H128.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H325.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H475.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H200.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H117.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H950.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H250.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H170.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H160.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H900.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H127.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H800.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H1000.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H500.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H135.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H375.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H190.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H119.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H129.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H180.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H750.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H275.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H118.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H650.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H121.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H120.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H122.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H220.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H400.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H850.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H115.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H525.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H150.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H600.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H123.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H575.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H425.root");  ggTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H140.root");  ggTree->Add(nameFile2);
    } 
  }
  
  TChain* VBFTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH125.root");    VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH125.root");    VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH125.root");    VBFTree->Add(nameFile2);
    if (moreFiles) {
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH115.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH210.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH230.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH115.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH210.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH230.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH115.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH210.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH230.root");  VBFTree->Add(nameFile2);
    }
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH125.root");  VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH125.root");  VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH125.root");  VBFTree->Add(nameFile2);
    if (moreFiles) {
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH128.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH135.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH850.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH525.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH750.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH126.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH116.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH117.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH121.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH145.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH127.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH122.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH129.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH119.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH550.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH118.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH124.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH123.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH128.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH135.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH850.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH525.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH750.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH126.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH116.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH117.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH121.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH145.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH127.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH122.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH129.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH119.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH550.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH118.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH124.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH123.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH200.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH150.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH128.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH135.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH700.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH950.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH650.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH500.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH300.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH600.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH450.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH375.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH250.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH850.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH275.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH525.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH750.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH350.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH126.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH160.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH116.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH220.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH117.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH575.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH190.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH425.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH325.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH475.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH121.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH145.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH127.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH122.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH129.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH130.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH119.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH1000.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH120.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH550.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH400.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH170.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH140.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH180.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH118.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH900.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH124.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH800.root");  VBFTree->Add(nameFile2);
      sprintf(nameFile2,"root://lxcms02//data//Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH123.root");  VBFTree->Add(nameFile2);
    }
  }

  TChain* VHTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    if (moreFiles) {
      
    }
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VH125.root");   VHTree->Add(nameFile2);
    if (moreFiles) {
      
    }
  }
  
  TChain* zzTree = new TChain(nameFile);
  if (also7TeV) {
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo4e.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo4mu.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2e2mu.root"); 
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo4tau.root");
  } else {
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4e.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4mu.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2mu.root"); 
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo4tau.root");
  }
  
  TChain* ggzzTree = new TChain(nameFile);
  if (also7TeV) {
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  } else {
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  }
  
  TChain* dataTree = new TChain(nameFile);
  if (also7TeV) {
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleEle.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleMu.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleOr.root");
  } else {
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleEle.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleMu.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleOr.root"); 
  }

  TChain* crTree = new TChain(nameFile);
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEMMssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMEEssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMMMssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMMMssTree.root"); 
  if (also7TeV) {
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEMMssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMEEssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMMMssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMMMssTree.root");
  }

  TChain* crosTree = new TChain(nameFile);
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEMMosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMEEosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMMMosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMMMosTree.root");	
  if (also7TeV) {
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEMMosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMEEosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMMMosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMMMosTree.root");	
  }

  TChain* crzjTree = new TChain(nameFile);
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEEEssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMEEssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMMMssTree.root");

  float mgg, mVBF, mzz, mggzz, mdata, mcr, mcros, mzj, mvh;
  float wgg, wVBF, wzz, wggzz, wzj, wvh;
  float ptgg, ptVBF, ptzz, ptggzz, ptdata, ptcr, ptzj, ptvh;
  int njgg, njVBF, njzz, njggzz, njdata, njcr, njzj, njvh;
  float nlogg = -999.;
  float nloVBF = -999.;
  float nlozz = -999.; 
  float nloggzz = -999.;
  float nlodata = -999.; 
  float nlocr = -999.; 
  float nlozj = -999.; 
  float nlovh = -999.; 
  float genptgg, genptvh, genptvbf;
  int gprIdvh; 

  ggTree->SetBranchAddress("ZZMass",&mgg);
  ggTree->SetBranchAddress("MC_weight",&wgg);
  ggTree->SetBranchAddress("ZZPt",&ptgg);
  ggTree->SetBranchAddress("NJets",&njgg);
  ggTree->SetBranchAddress("GenHPt",&genptgg);
  if (withNLOMela) ggTree->SetBranchAddress("melaLDWithPtY",&nlogg);
  else ggTree->SetBranchAddress("ZZLD",&nlogg);

  VBFTree->SetBranchAddress("ZZMass",&mVBF);
  VBFTree->SetBranchAddress("MC_weight",&wVBF);
  VBFTree->SetBranchAddress("ZZPt",&ptVBF);
  VBFTree->SetBranchAddress("NJets",&njVBF);
  VBFTree->SetBranchAddress("GenHPt",&genptvbf);
  if (withNLOMela) VBFTree->SetBranchAddress("melaLDWithPtY",&nloVBF);
  else VBFTree->SetBranchAddress("ZZLD",&nloVBF);  

  VHTree->SetBranchAddress("ZZMass",&mvh);
  VHTree->SetBranchAddress("MC_weight",&wvh);
  VHTree->SetBranchAddress("ZZPt",&ptvh);
  VHTree->SetBranchAddress("NJets",&njvh);
  VHTree->SetBranchAddress("genProcessId",&gprIdvh);
  VHTree->SetBranchAddress("GenHPt",&genptvh);
  if (withNLOMela) VHTree->SetBranchAddress("melaLDWithPtY",&nlovh);
  else VHTree->SetBranchAddress("ZZLD",&nlovh);

  zzTree->SetBranchAddress("ZZMass",&mzz);
  zzTree->SetBranchAddress("MC_weight",&wzz);
  zzTree->SetBranchAddress("ZZPt",&ptzz);
  zzTree->SetBranchAddress("NJets",&njzz);
  if (withNLOMela) zzTree->SetBranchAddress("melaLDWithPtY",&nlozz);
  else zzTree->SetBranchAddress("ZZLD",&nlozz);

  ggzzTree->SetBranchAddress("ZZMass",&mggzz);
  ggzzTree->SetBranchAddress("MC_weight",&wggzz);
  ggzzTree->SetBranchAddress("ZZPt",&ptggzz);
  ggzzTree->SetBranchAddress("NJets",&njggzz);
  if (withNLOMela) ggzzTree->SetBranchAddress("melaLDWithPtY",&nloggzz);
  else ggzzTree->SetBranchAddress("ZZLD",&nloggzz);

  dataTree->SetBranchAddress("ZZMass",&mdata);
  dataTree->SetBranchAddress("ZZPt",&ptdata);
  dataTree->SetBranchAddress("NJets",&njdata);
  if (withNLOMela) dataTree->SetBranchAddress("melaLDWithPtY",&nlodata);
  else dataTree->SetBranchAddress("ZZLD",&nlodata);

  crTree->SetBranchAddress("ZZMass",&mcr);
  crTree->SetBranchAddress("ZZPt",&ptcr);
  crTree->SetBranchAddress("NJets",&njcr);
  if (withNLOMela) crTree->SetBranchAddress("melaLDWithPtY",&nlocr);
  else crTree->SetBranchAddress("ZZLD",&nlocr);

  crosTree->SetBranchAddress("ZZMass",&mcros);

  crzjTree->SetBranchAddress("ZZMass",&mzj);
  crzjTree->SetBranchAddress("MC_weight",&wzj);
  crzjTree->SetBranchAddress("ZZPt",&ptzj);
  crzjTree->SetBranchAddress("NJets",&njzj);
  if (withNLOMela) crzjTree->SetBranchAddress("melaLDWithPtY",&nlozj);
  else crzjTree->SetBranchAddress("ZZLD",&nlozj);

  char nameSyst[200]; 
  char UcasePt[8] = "PT";
  if (overM == 1) sprintf(UcasePt,"PToverM");
  else if (overM == 2) sprintf(UcasePt,"lnPToverM");
  char LcasePt[8] = "pt";
  if (overM == 1) sprintf(LcasePt,"ptoverm");
  else if (overM == 2) sprintf(LcasePt,"lnptoverm");

  TH1F* wHalf[nWeightResum];
  TH1F* wQuar[nWeightResum];
  TH1F* wOne[nWeightResum];

  TFile f125u("newweights/weightResumUp_125GeV_8TeV.root");
  wOne[0] = (TH1F*)f125u.Get("wH");
  TFile f125d("newweights/weightResumDown_125GeV_8TeV.root");
  wQuar[0] = (TH1F*)f125d.Get("wH");
  TFile f400("newweights/weightResumDef_400GeV_8TeV.root");
  wHalf[1] = (TH1F*)f400.Get("wH");
  TFile f400u("newweights/weightResumUp_400GeV_8TeV.root");
  wOne[1] = (TH1F*)f400u.Get("wH");
  TFile f400d("newweights/weightResumDown_400GeV_8TeV.root");
  wQuar[1] = (TH1F*)f400d.Get("wH");
  
  sprintf(nameFile,"newweights/ggH125_POWHEG_btMassErrors.root");
  TFile pwhgNoMass(nameFile);
  TH1F* wbtMass = (TH1F*)pwhgNoMass.Get("wei");

  TH1F* wNLO[nWeightNLOvh];

  TFile f115w("newweights/wh115_weightsNLO.root");
  wNLO[0] = (TH1F*)f115w.Get("wei");
  TFile f125w("newweights/wh125_weightsNLO.root");
  wNLO[1] = (TH1F*)f125w.Get("wei");
  TFile f140w("newweights/wh140_weightsNLO.root");
  wNLO[2] = (TH1F*)f140w.Get("wei");

  // cout << wNLO[1]->GetEntries() << endl;

  TCanvas can2("can2","The 2nd canvas",15.,15.,700.,500.); 
  TCanvas can("can","The canvas",5.,5.,600.,700.); 
  can.Divide(1,2);

  float massLimits[massRanges+1] = {100.,150.,200.,600.,1600.};   // fill now  
  float melaLimits[melaRanges+1] = {0.0,0.3,0.6,1.0}; 
  
  // gg 
  TH2F* pth = new TH2F("pth","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  pth->Sumw2();

  TH2F* pthMass[massRanges];
  TH2F* pthMela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"pthMass%d",i);
    pthMass[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    pthMass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"pthMela%d",i);
    pthMela[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    pthMela[i]->Sumw2();
  }

  // vbf
  TH2F* ptvbf = new TH2F("ptvbf","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  ptvbf->Sumw2();
  
  TH2F* ptvbfMass[massRanges];
  TH2F* ptvbfMela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"ptvbfMass%d",i);
    ptvbfMass[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    ptvbfMass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"ptvbfMela%d",i);
    ptvbfMela[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    ptvbfMela[i]->Sumw2();
  }

  // zz
  TH2F* pth1 = new TH2F("pth1","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  pth1->Sumw2();
  TH2F* pth1Mass[massRanges];
  TH2F* pth1Mela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"pth1Mass%d",i);
    pth1Mass[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    pth1Mass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"pth1Mela%d",i);
    pth1Mela[i] = new TH2F(nameFile,"pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
    pth1Mela[i]->Sumw2();
  }

  // other
  TH2F* pth2 = new TH2F("pth2","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  TH2F* pth3 = new TH2F("pth3","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  TH2F* pth4 = new TH2F("pth4","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  TH2F* pth5 = new TH2F("pth5","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  TH2F* pth6 = new TH2F("pth6","pt",mbins,mbinLimits,ptmbins,ptmbinLimits);
  pth2->Sumw2();
  pth3->Sumw2();
  pth4->Sumw2();
  pth5->Sumw2();
  pth6->Sumw2();

  // for single Z data
  const int nbins = 10;
  float binlimits[nbins+1] = {0.,10.,20.,30.,40.,50.,70.,90.,110.,150.,190.}; 
  TH1F* ptalt = new TH1F("ptalt","pt",nbins,binlimits);
  TH1F* ptalt1 = new TH1F("ptalt1","pt",nbins,binlimits);
  ptalt->Sumw2();
  ptalt1->Sumw2();
 
  if (whichtype == -7) {

    for (Int_t iEvt = 0; iEvt < VHTree->GetEntries() ; ++iEvt) {
      VHTree->GetEntry(iEvt);
      if (mvh < massLimits[massRanges] && mvh > massLimits[0] && ptvh < 400. && notVBFtagged(njvh)) {
	if (gprIdvh == 26) {
          int theHist = returnClosestMass(nWeightNLOvh,massNLOvh,mvh);   
	  int theBin = wNLO[theHist]->FindBin(genptvh);
	  pth->Fill(mvh,ptVar(ptvh,mvh,overM),wvh*wNLO[theHist]->GetBinContent(theBin));
	  pth2->Fill(mvh,ptVar(ptvh,mvh,overM),wvh);
	  ptvbf->Fill(mvh,ptVar(ptvh,mvh,overM),wvh);
	}
      }
    } 
    
    sprintf(nameSyst,"NLOLO_ZH");
    
  }
  
  if (whichtype == -6) {
    
    for (Int_t iEvt = 0; iEvt < VHTree->GetEntries() ; ++iEvt) {
      VHTree->GetEntry(iEvt);
      if (mvh < massLimits[massRanges] && mvh > massLimits[0] && ptvh < 400. && notVBFtagged(njvh)) {
	if (gprIdvh == 24) {
          int theHist = returnClosestMass(nWeightNLOvh,massNLOvh,mvh);   
	  int theBin = wNLO[theHist]->FindBin(genptvh);
	  pth->Fill(mvh,ptVar(ptvh,mvh,overM),wvh*wNLO[theHist]->GetBinContent(theBin));
	  pth2->Fill(mvh,ptVar(ptvh,mvh,overM),wvh);
	  ptvbf->Fill(mvh,ptVar(ptvh,mvh,overM),wvh);
	}
      }
    } 
    
    sprintf(nameSyst,"NLOLO_WH");
    
  }

  if (whichtype == -5) {

    TH1F* pdfh[3];
    TFile* pdfsfile0 = TFile::Open("newweights/VBF_standard.root");
    sprintf(nameFile,"pdfh0");
    pdfh[0] = (TH1F*)((TH2F*)pdfsfile0->Get("Pt_sig"))->ProjectionY(nameFile);
    pdfh[0]->Rebin(4);   pdfh[0]->Sumw2();
    pdfh[0]->GetXaxis()->SetRange(1,80);
    pdfh[0]->Scale(1./pdfh[0]->Integral());

    TFile* pdfsfile[3];
    pdfsfile[1] = TFile::Open("newweights/VBF_MSTW.root");
    pdfsfile[2] = TFile::Open("newweights/VBF_NNPDF.root");
   
    TH1F* ratiopdf[2];
    for (int i = 1; i < 3; i++) {
      sprintf(nameFile,"pdfh%d",i);
      pdfh[i] = (TH1F*)((TH2F*)pdfsfile[i]->Get("Pt_sig"))->ProjectionY(nameFile);
      pdfh[i]->Rebin(4);   pdfh[i]->Sumw2();
      pdfh[i]->GetXaxis()->SetRange(1,80);
      pdfh[i]->Scale(1./pdfh[i]->Integral());
      ratiopdf[i-1] = (TH1F*)pdfh[i]->Clone();
      ratiopdf[i-1]->Divide(pdfh[i],pdfh[0]);
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	int theBin = pdfh[0]->FindBin(genptvbf);
 	pth->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
	pth1->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*ratiopdf[0]->GetBinContent(theBin));
      	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*thisDiff);
        // cout << "thisDiff PDF-VBF = " << thisDiff << endl;   
      }
    } 

    sprintf(nameSyst,"PDF-VBF");
  }

  if (whichtype == -4) {

    TH1F* scaleh[7];
    TFile* scalesfile0 = TFile::Open("newweights/VBF_standard.root");
    sprintf(nameFile,"scaleh0");
    scaleh[0] = (TH1F*)((TH2F*)scalesfile0->Get("Pt_sig"))->ProjectionY(nameFile);
    scaleh[0]->Rebin(4);   scaleh[0]->Sumw2();
    scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[0]->Scale(1./scaleh[0]->Integral());

    TFile* scalesfile[7];
    scalesfile[1] = TFile::Open("newweights/VBF_mufdouble.root");
    scalesfile[2] = TFile::Open("newweights/VBF_mufhalf.root");
    scalesfile[3] = TFile::Open("newweights/VBF_murdouble.root");
    scalesfile[4] = TFile::Open("newweights/VBF_murhalf.root");
    scalesfile[5] = TFile::Open("newweights/VBF_mufdouble_murdouble.root");
    scalesfile[6] = TFile::Open("newweights/VBF_mufhalf_murhalf.root");
    
    TH1F* ratioscale[6];
    for (int i = 1; i < 7; i++) {
      sprintf(nameFile,"scaleh%d",i);
      scaleh[i] = (TH1F*)((TH2F*)scalesfile[i]->Get("Pt_sig"))->ProjectionY(nameFile);
      scaleh[i]->Rebin(4);   scaleh[i]->Sumw2();
      scaleh[i]->GetXaxis()->SetRange(1,80);
      scaleh[i]->Scale(1./scaleh[i]->Integral());
      ratioscale[i-1] = (TH1F*)scaleh[i]->Clone();
      ratioscale[i-1]->Divide(scaleh[i],scaleh[0]); 
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	int theBin = scaleh[0]->FindBin(genptvbf);
 	pth->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
	pth1->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*ratioscale[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*thisDiff); 
        // cout << "thisDiff scale-VBF = " << thisDiff << endl; 
      }
    } 

    sprintf(nameSyst,"scale-VBF");
    
  }  

  if (whichtype == -3) {

    sprintf(nameSyst,"TopMass");
    
    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
        float oldW = 1.;
        if (mgg > 400) {
	  int theHist = returnClosestMass(nWeightResum,massResum,mgg);
	  oldW = wHalf[theHist]->GetBinContent(wHalf[theHist]->FindBin(genptgg));
	}
	float newW = 1.;
      	newW = wbtMass->GetBinContent(wbtMass->FindBin(genptgg));
	pth->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
	pth2->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*newW);
	ptvbf->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*newW);
      }
    } 
    
  } 
  
  if (whichtype == -2) {
    
    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
        int theHist = returnClosestMass(nWeightResum,massResum,mgg); 
	int theBin = wQuar[theHist]->FindBin(genptgg);
	float oldW = 1.;
        if (mgg > 400) oldW = wHalf[theHist]->GetBinContent(theBin);
	pth->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
	pth1->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*wQuar[theHist]->GetBinContent(theBin));
	pth2->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*wOne[theHist]->GetBinContent(theBin));
	ptvbf->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*wQuar[theHist]->GetBinContent(theBin));
      }
    } 
    
    sprintf(nameSyst,"Resummation");
    
  }
  
  if (whichtype == -1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	float oldW = 1.;
        if (mgg > 400) {
	  int theHist = returnClosestMass(nWeightResum,massResum,mgg);
	  oldW = wHalf[theHist]->GetBinContent(wHalf[theHist]->FindBin(genptgg));
	}
	pth->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
	pth1->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*0.84);  // gg: +/-16%
	pth2->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*1.16);
	ptvbf->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW*1.16);  
      }
    } 
    
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
        pth1->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*1.02);  // VBF: +/-2%
	pth2->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*0.98);
	ptvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF*0.98);
	
      }
    }
    
    sprintf(nameSyst,"VBFfraction");
    
  }
  
  if (whichtype == 0) {
    
    if (overM == 0) sprintf(nameFile2,"p_{T} [GeV/c]");
    else if (overM == 1) sprintf(nameFile2,"p_{T}/m_{4l}");
    else sprintf(nameFile2,"ln(p_{T}/m_{4l})");

    TH2F* ptkdgg = new TH2F("ptkdgg","pt vs mela",16,0.,1.,ptmbins,ptmbinLimits);
    ptkdgg->Sumw2();
    TH2F* ptmhgg = new TH2F("ptmhgg","pt vs m",16,100.,180.,ptmbins,ptmbinLimits);
    ptmhgg->Sumw2();
    
    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      float oldW = 1.;
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
        if (mgg > 400) {
	  int theHist = returnClosestMass(nWeightResum,massResum,mgg);
	  oldW = wHalf[theHist]->GetBinContent(wHalf[theHist]->FindBin(genptgg));
	}
	pth->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
        ptkdgg->Fill(nlogg,ptVar(ptgg,mgg,overM),wgg*oldW); 
	for (int i = 0; i < melaRanges; i++) {
	  if (nlogg > melaLimits[i] && nlogg < melaLimits[i+1]) 
	    pthMela[i]->Fill(mgg,ptVar(ptgg,mgg,overM));
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mgg > massLimits[i] && mgg < massLimits[i+1] && ptgg < 400. && notVBFtagged(njgg)) 
	  pthMass[i]->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
      }      
      ptmhgg->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
    } 
    
    sprintf(nameFile,"selRootFiles/%s_gg_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_gg_TEMPL_7TeV.root",UcasePt);
    TFile f(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth->SetName(nameFile);
    // pth->Scale(1./pth->Integral());
    adjustHistogram(pth);   pth->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      pthMass[i]->SetName(nameFile);
      // pthMass[i]->Scale(1./pthMass[i]->Integral());
      adjustHistogram(pthMass[i]);   pthMass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      pthMela[i]->SetName(nameFile);
      // pthMela[i]->Scale(1./pthMela[i]->Integral());
      adjustHistogram(pthMela[i]);   pthMela[i]->Write();
    }
    
    ptmhgg->SetName("ptVsM");   ptmhgg->Write();
    TProfile* tpgg = (TProfile*)ptmhgg->ProfileX();
    tpgg->SetName("ptVsMProf");   tpgg->Write();

    ptkdgg->SetName("ptVsMELA");   ptkdgg->Write();
    TProfile* tpkdgg = (TProfile*)ptkdgg->ProfileX();
    tpkdgg->SetName("ptVsMELAProf");   tpkdgg->Write();

    // evalBinMigration(pth,pth,"",false);   //reinsert dop
    f.Close();
    can2.cd();   gPad->SetLogx();
    pth->GetYaxis()->SetTitle(nameFile2);
    pth->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_ggDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_ggDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    TH2F* ptkdvbf = new TH2F("ptkdvbf","pt vs mela",16,0.,1.,ptmbins,ptmbinLimits);
    ptkdvbf->Sumw2();
    TH2F* ptmhvbf = new TH2F("ptmhvbf","pt vs m",16,100.,180.,ptmbins,ptmbinLimits);
    ptmhvbf->Sumw2(); 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	ptvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
        ptkdvbf->Fill(nloVBF,ptVar(ptVBF,mVBF,overM),wVBF);
	for (int i = 0; i < melaRanges; i++) {
	if (nloVBF > melaLimits[i] && nloVBF < melaLimits[i+1]) 
	  ptvbfMela[i]->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mVBF > massLimits[i] && mVBF < massLimits[i+1] && ptVBF < 400. && notVBFtagged(njVBF)) 
	  ptvbfMass[i]->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
      }
      ptmhvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
    }

    sprintf(nameFile,"selRootFiles/%s_vbf_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_vbf_TEMPL_7TeV.root",UcasePt);
    TFile f1(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    ptvbf->SetName(nameFile);
    // ptvbf->Scale(1./ptvbf->Integral());
    adjustHistogram(ptvbf);   ptvbf->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      ptvbfMass[i]->SetName(nameFile);
      // ptvbfMass[i]->Scale(1./ptvbfMass[i]->Integral());
      adjustHistogram(ptvbfMass[i]);   ptvbfMass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      ptvbfMela[i]->SetName(nameFile);
      // ptvbfMela[i]->Scale(1./ptvbfMela[i]->Integral());
      adjustHistogram(ptvbfMela[i]);   ptvbfMela[i]->Write();
    }
    
    ptmhvbf->SetName("ptVsM");   ptmhvbf->Write();
    TProfile* tpvbf = (TProfile*)ptmhvbf->ProfileX();
    tpvbf->SetName("ptVsMProf");   tpvbf->Write();

    ptkdvbf->SetName("ptVsMELA");   ptkdvbf->Write();
    TProfile* tpkdvbf = (TProfile*)ptkdvbf->ProfileX();
    tpkdvbf->SetName("ptVsMELAProf");   tpkdvbf->Write();

    // evalBinMigration(ptvbf,ptvbf,"",false);
    f1.Close();
    can2.cd();   gPad->SetLogx();
    ptvbf->GetYaxis()->SetTitle(nameFile2);
    ptvbf->GetXaxis()->SetTitle("m_{4l} [GeV]");
    ptvbf->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_vbfDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_vbfDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    TH2F* ptkdzz = new TH2F("ptkdzz","pt vs mela",16,0.,1.,ptmbins,ptmbinLimits);
    ptkdzz->Sumw2();
    TH2F* ptmhzz = new TH2F("ptmhzz","pt vs m",16,100.,180.,ptmbins,ptmbinLimits);
    ptmhzz->Sumw2();

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth1->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
        ptkdzz->Fill(nlozz,ptVar(ptzz,mzz,overM),wzz);
	for (int i = 0; i < melaRanges; i++) {
	  if (nlozz > melaLimits[i] && nlozz < melaLimits[i+1]) 
	    pth1Mela[i]->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mzz > massLimits[i]-7. && mzz < massLimits[i+1]+7. && ptzz < 400. && notVBFtagged(njzz)) 
	  pth1Mass[i]->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
      }
      ptmhzz->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
    } 

    sprintf(nameFile,"selRootFiles/%s_zz_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_zz_TEMPL_7TeV.root",UcasePt);
    TFile f2(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth1->SetName(nameFile);
    // pth1->Scale(1./pth1->Integral());
    adjustHistogram(pth1);   pth1->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      pth1Mass[i]->SetName(nameFile);
      // pth1Mass[i]->Scale(1./pth1Mass[i]->Integral());
      adjustHistogram(pth1Mass[i]);   pth1Mass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      pth1Mela[i]->SetName(nameFile);
      // pth1Mela[i]->Scale(1./pth1Mela[i]->Integral());
      adjustHistogram(pth1Mela[i]);   pth1Mela[i]->Write();
    }
    
    ptmhzz->SetName("ptVsM");   ptmhzz->Write();
    TProfile* tpzz = (TProfile*)ptmhzz->ProfileX();
    tpzz->SetName("ptVsMProf");   tpzz->Write();

    ptkdzz->SetName("ptVsMELA");   ptkdzz->Write();
    TProfile* tpkdzz = (TProfile*)ptkdzz->ProfileX();
    tpkdzz->SetName("ptVsMELAProf");   tpkdzz->Write();

    //    evalBinMigration(pth1,pth1,"",false);
    f2.Close();
    can2.cd();   gPad->SetLogx();
    pth1->GetYaxis()->SetTitle(nameFile2);
    pth1->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth1->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_zzDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_zzDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    TH2F* ptkdcr = new TH2F("ptkdcr","pt vs mela",16,0.,1.,ptmbins,ptmbinLimits);
    ptkdcr->Sumw2();
    TH2F* ptmhcr = new TH2F("ptmhcr","pt vs m",16,100.,180.,ptmbins,ptmbinLimits);
    ptmhcr->Sumw2();

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < massLimits[massRanges] && mcr > massLimits[0] && ptcr < 400. && notVBFtagged(njcr)) {
	pth2->Fill(mcr,ptVar(ptcr,mcr,overM));
        ptkdcr->Fill(nlocr,ptVar(ptcr,mcr,overM));
      }
      ptmhcr->Fill(mcr,ptVar(ptcr,mcr,overM));
    } 

    sprintf(nameFile,"selRootFiles/%s_zx_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_zx_TEMPL_7TeV.root",UcasePt);
    TFile f3(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth2->SetName(nameFile);
    // pth2->Scale(1./pth2->Integral());
    adjustHistogram(pth2);   pth2->Write();

    ptmhcr->SetName("ptVsM");   ptmhcr->Write();
    TProfile* tpcr = (TProfile*)ptmhcr->ProfileX();
    tpcr->SetName("ptVsMProf");   tpcr->Write();

    ptkdcr->SetName("ptVsMELA");   ptkdcr->Write();
    TProfile* tpkdcr = (TProfile*)ptkdcr->ProfileX();
    tpkdcr->SetName("ptVsMELAProf");   tpkdcr->Write();

    //    evalBinMigration(pth2,pth2,"",false);
    f3.Close();
    can2.cd();   gPad->SetLogx();
    pth2->GetYaxis()->SetTitle(nameFile2);
    pth2->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth2->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_zxDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_zxDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    for (Int_t iEvt = 0; iEvt < ggzzTree->GetEntries() ; ++iEvt) {
      ggzzTree->GetEntry(iEvt);
      if (mggzz < massLimits[massRanges] && mggzz > massLimits[0] && ptggzz < 400. && notVBFtagged(njggzz)) {
	pth3->Fill(mggzz,ptVar(ptggzz,mggzz,overM),wggzz);
      }
    } 

    sprintf(nameFile,"selRootFiles/%s_ggzz_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_ggzz_TEMPL_7TeV.root",UcasePt);
    TFile f4(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth3->SetName(nameFile);
    // pth3->Scale(1./pth3->Integral());
    adjustHistogram(pth3);   pth3->Write();
    // evalBinMigration(pth3,pth3,"",false);
    f4.Close();
    can2.cd();   gPad->SetLogx();
    pth3->GetYaxis()->SetTitle(nameFile2);
    pth3->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth3->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_ggzzDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_ggzzDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    for (Int_t iEvt = 0; iEvt < VHTree->GetEntries() ; ++iEvt) {
      VHTree->GetEntry(iEvt);
      if (mvh < massLimits[massRanges] && mvh > massLimits[0] && ptvh < 400. && notVBFtagged(njvh)) {
        int theHist = returnClosestMass(nWeightNLOvh,massNLOvh,mvh);  
      	int theBin = wNLO[theHist]->FindBin(genptvh);
      	if (gprIdvh == 24) {
	  pth4->Fill(mvh,ptVar(ptvh,mvh,overM),wvh*wNLO[theHist]->GetBinContent(theBin));
	} else if (gprIdvh == 26) {
	  pth5->Fill(mvh,ptVar(ptvh,mvh,overM),wvh*wNLO[theHist]->GetBinContent(theBin));
	} else {
	  pth6->Fill(mvh,ptVar(ptvh,mvh,overM),wvh);
	}
      }
    } 
    
    sprintf(nameFile,"selRootFiles/%s_wh_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_wh_TEMPL_7TeV.root",UcasePt);
    TFile f5(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth4->SetName(nameFile);
    // pth4->Scale(1./pth4->Integral());
    adjustHistogram(pth4);   pth4->Write();
    // evalBinMigration(pth4,pth4,"",false);
    f5.Close();
    can2.cd();   gPad->SetLogx();
    pth4->GetYaxis()->SetTitle(nameFile2);
    pth4->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth4->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_whDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_whDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    sprintf(nameFile,"selRootFiles/%s_zh_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_zh_TEMPL_7TeV.root",UcasePt);
    TFile f6(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth5->SetName(nameFile);
    // pth5->Scale(1./pth5->Integral());
    adjustHistogram(pth5);   pth5->Write();
    // evalBinMigration(pth5,pth5,"",false);
    f6.Close();
    can2.cd();   gPad->SetLogx();
    pth5->GetYaxis()->SetTitle(nameFile2);
    pth5->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth5->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_zhDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_zhDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

    sprintf(nameFile,"selRootFiles/%s_tth_TEMPL_8TeV.root",UcasePt);
    if (also7TeV) sprintf(nameFile,"selRootFiles/%s_tth_TEMPL_7TeV.root",UcasePt);
    TFile f7(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth6->SetName(nameFile);
    // pth6->Scale(1./pth6->Integral());
    adjustHistogram(pth6);   pth6->Write();
    // evalBinMigration(pth6,pth6,"",false);
    f7.Close();
    can2.cd();   gPad->SetLogx();
    pth6->GetYaxis()->SetTitle(nameFile2);
    pth6->GetXaxis()->SetTitle("m_{4l} [GeV]");
    pth6->Draw("COLZ");
    sprintf(nameFile,"newfigs/%s_tthDefaultTemplate_%dTeV.gif",LcasePt,sqrts);
    can2.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_tthDefaultTemplate_%dTeV.pdf",LcasePt,sqrts);
    can2.SaveAs(nameFile); 

  }

  
  if (whichtype == 1) {
    
    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
        float oldW = 1.;
        if (mgg > 400) {
	  int theHist = returnClosestMass(nWeightResum,massResum,mgg);
	  oldW = wHalf[theHist]->GetBinContent(wHalf[theHist]->FindBin(genptgg));
	}
	pth1->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*oldW);
      }
    } 
    
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth1->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
      }
    }
    
    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth2->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
      }
    } 
    
    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < massLimits[massRanges] && mcr > massLimits[0] && ptcr < 400. && notVBFtagged(njcr)) {
	pth->Fill(mcr,ptVar(ptcr,mcr,overM));
      }
    } 
    
    sprintf(nameSyst,"ZplusX");
    
  }

  if (whichtype == 2) {
    
    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (ptcr < 400. && notVBFtagged(njcr)) {          // do it in the whole mass range
	pth->Fill(mcr,ptVar(ptcr,mcr,overM));
      }
    } 
    
    for (Int_t iEvt = 0; iEvt < crzjTree->GetEntries() ; ++iEvt) {
      crzjTree->GetEntry(iEvt);
      if (ptzj < 400. && notVBFtagged(njzj)) {
	pth2->Fill(mzj,ptVar(ptzj,mzj,overM),wzj);
	ptvbf->Fill(mzj,ptVar(ptzj,mzj,overM),wzj);
      }
    }
    
    sprintf(nameSyst,"ZjetsXcheck");
    
  }

  if (whichtype == 3) {
    
    float unbMaxBin = 200.;
    if (overM) unbMaxBin = 0.8;
    
    TH1F* dataUnb = new TH1F("dataUnb","data",20,0.,unbMaxBin);
    TH1F* zzUnb = new TH1F("zzUnb","zz",20,0.,unbMaxBin);
    TH1F* zxUnb = new TH1F("zxUnb","zx",20,0.,unbMaxBin);
    dataUnb->Sumw2();
    zzUnb->Sumw2();
    zxUnb->Sumw2();

    for (Int_t iEvt = 0; iEvt < dataTree->GetEntries() ; ++iEvt) {
      dataTree->GetEntry(iEvt);
      if (mdata < 400. && mdata > massLimits[massRanges] && ptdata < 200. && notVBFtagged(njdata)) {
	dataUnb->Fill(ptVar(ptdata,mdata,overM));
      }
    } 

    float lumi = 12.2;
    if (also7TeV) lumi = 5.1;

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 400. && mzz > massLimits[massRanges] && ptzz < 200. && notVBFtagged(njzz)) {
	zzUnb->Fill(ptVar(ptzz,mzz,overM),wzz*lumi);
      }
    } 

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < 400. && mcr > massLimits[massRanges] && ptcr < 200. && notVBFtagged(njcr)) {
	zxUnb->Fill(ptVar(ptcr,mcr,overM));
      }
    } 

    // Subtract fakes 
    float weightZX = (4.193/20.516)*(zzUnb->Integral()/zxUnb->Integral());
    dataUnb->Add(dataUnb,zxUnb,1,-weightZX);
    dataUnb->SetMarkerColor(1);    dataUnb->SetLineColor(1);    dataUnb->SetLineWidth(2);
    zzUnb->SetMarkerColor(2);   zzUnb->SetLineColor(2);   zzUnb->SetLineWidth(2);
    // zxUnb->SetMarkerColor(3);   zxUnb->SetLineColor(3);   zxUnb->SetLineWidth(2);

    can.cd(1);
    gPad->SetBottomMargin(0.0);
    dataUnb->GetXaxis()->SetLabelColor(kWhite);
    dataUnb->GetYaxis()->SetTitle("Entries");
    dataUnb->Draw();
    zzUnb->Draw("SAME");
    // zxUnb->Draw("SAME");

    TH1F* diffpt = (TH1F*)zxUnb->Clone();
    TH1F* pullpt = (TH1F*)zxUnb->Clone();
    diffpt->Add(dataUnb,zzUnb,1,-1);
    pullpt->Divide(diffpt,dataUnb);
    // remove outliers
    float xp[20];     float yp[20];
    float xperr[20];  float yperr[20];
    for (Int_t iBin = 1; iBin <= pullpt->GetNbinsX() ; ++iBin) {
      xperr[iBin-1] = 0.;
      if (fabs(pullpt->GetBinContent(iBin)) > 2.) {
	xp[iBin-1] = 999.;   yp[iBin-1] = 0.;     yperr[iBin-1] = 0.;
      } else {
	xp[iBin-1] = pullpt->GetBinCenter(iBin);   
	yp[iBin-1] = pullpt->GetBinContent(iBin);     
	yperr[iBin-1] = pullpt->GetBinError(iBin);;
      }
    }
    TGraphErrors* pullg = new TGraphErrors(20,xp,yp,xperr,yperr);

    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-2.);
    pullpt->SetMaximum(2.);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
    TF1 *mypol1 = new TF1("mypol1","pol1");
    pullg->Fit("mypol1","","",0.,9.*unbMaxBin/20.);
    pullg->Draw("PSAME");

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*(1. + mypol1->Eval(ptVar(ptzz,mzz,overM))));
      }
    } 
    
    sprintf(nameSyst,"UnbRegion");

    sprintf(nameFile,"newfigs/%s_systUnbRegion.gif",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_systUnbRegion.pdf",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    // return;

  }

  if (whichtype == 4) {
    
    float dataSz[nbins] = {465.0,251.5,116.0,59.8,33.8,18.1,7.79,4.75,1.93,0.60};
    float dataErr[nbins] = {7.5,5.4,4.0,2.7,1.8,0.9,0.54,0.42,0.17,0.10};
    float powSz[nbins] = {454.3,256.7,118.0,58.0,32.0,18.2,7.6,3.4,1.2,0.59};
    float powErr[nbins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    for (Int_t iBin = 0; iBin < nbins; ++iBin) {
      ptalt->SetBinContent(iBin+1,dataSz[iBin]);
      ptalt->SetBinError(iBin+1,dataErr[iBin]);
      ptalt1->SetBinContent(iBin+1,powSz[iBin]);
      ptalt1->SetBinError(iBin+1,powErr[iBin]);
    }

    ptalt->SetMarkerColor(1);    ptalt->SetLineColor(1);    ptalt->SetLineWidth(2);
    ptalt1->SetMarkerColor(2);   ptalt1->SetLineColor(2);   ptalt1->SetLineWidth(2);
    can.cd(1);
    gPad->SetBottomMargin(0.0);
    ptalt->GetXaxis()->SetLabelColor(kWhite);
    ptalt->GetYaxis()->SetTitle("Entries");
    ptalt->Draw();
    ptalt1->Draw("SAME");

    TH1F* diffpt = (TH1F*)ptalt->Clone();
    TH1F* pullpt = (TH1F*)ptalt->Clone();
    diffpt->Add(ptalt,ptalt1,1,-1);
    pullpt->Divide(diffpt,ptalt);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-1.);
    pullpt->SetMaximum(1.);
    pullpt->GetXaxis()->SetLabelColor(kBlack);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
    TF1 *mypol1 = new TF1("mypol1","pol1");
    pullpt->Fit("mypol1","","");

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*(1. + mypol1->Eval(ptzz)));
        pth2->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
        ptvbf->Fill(mzz,ptVar(ptzz,mzz,overM),wzz); 
      }
    } 

    sprintf(nameSyst,"SingleZ");    
    
    sprintf(nameFile,"newfigs/%s_systSingleZ.gif",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_systSingleZ.pdf",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);

  }

  if (whichtype == 5) {

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
   	pth2->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz);
      }
    } 

    for (Int_t iEvt = 0; iEvt < ggzzTree->GetEntries() ; ++iEvt) {
      ggzzTree->GetEntry(iEvt);
      if (mggzz < massLimits[massRanges] && mggzz > massLimits[0] && ptggzz < 400. && notVBFtagged(njggzz)) {
	pth1->Fill(mggzz,ptVar(ptggzz,mggzz,overM),wggzz);
 	pth->Fill(mggzz,ptVar(ptggzz,mggzz,overM),wggzz);
       }
    } 

    sprintf(nameSyst,"ggZZ");

  }

  if (whichtype == 6) {

    TH2F* pdfh[3];
    TH2F* ratiopdf[2];

    TFile* pdfsfile0 = TFile::Open("newweights/ZZ_standard.root");
    pdfh[0] = (TH2F*)pdfsfile0->Get("Pt_bkg");
    pdfh[0]->Rebin2D(4,4);   pdfh[0]->Sumw2();
    // pdfh[0]->GetXaxis()->SetRange(1,80);
    pdfh[0]->Scale(1./pdfh[0]->Integral());

    TFile* pdfsfile1 = TFile::Open("newweights/ZZ_MSTW.root");
    pdfh[1] = (TH2F*)pdfsfile1->Get("Pt_bkg");
    pdfh[1]->Rebin2D(4,4);   pdfh[1]->Sumw2();
    // pdfh[0]->GetXaxis()->SetRange(1,80);
    pdfh[1]->Scale(1./pdfh[1]->Integral());
    ratiopdf[0] = (TH2F*)pdfh[1]->Clone();
    ratiopdf[0]->Divide(pdfh[1],pdfh[0]);

    TFile* pdfsfile2 = TFile::Open("newweights/ZZ_NNPDF.root");
    pdfh[2] = (TH2F*)pdfsfile2->Get("Pt_bkg");
    pdfh[2]->Rebin2D(4,4);   pdfh[2]->Sumw2();
    // pdfh[0]->GetXaxis()->SetRange(1,80);
    pdfh[2]->Scale(1./pdfh[2]->Integral());
    ratiopdf[1] = (TH2F*)pdfh[2]->Clone();
    ratiopdf[1]->Divide(pdfh[2],pdfh[0]);         

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	int theBin = pdfh[0]->FindBin(mzz,ptzz);
 	pth->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	pth1->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*ratiopdf[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
        // cout << "thisDiff PDF-ZZ = " << thisDiff << endl; 
	ptvbf->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*thisDiff); 
      }
    } 
    
    sprintf(nameSyst,"PDF-ZZ");
  }
  
  if (whichtype == 7) {

    TH2F* scaleh[8];
    TH2F* ratioscale[7];

    TFile* scalesfile0 = TFile::Open("newweights/ZZ_standard.root");
    scaleh[0] = (TH2F*)scalesfile0->Get("Pt_bkg");
    scaleh[0]->Rebin2D(4,4);   scaleh[0]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[0]->Scale(1./scaleh[0]->Integral());

    TFile* scalesfile1 = TFile::Open("newweights/ZZ_mufdouble.root");
    scaleh[1] = (TH2F*)scalesfile1->Get("Pt_bkg");
    scaleh[1]->Rebin2D(4,4);   scaleh[1]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[1]->Scale(1./scaleh[1]->Integral());
    ratioscale[0] = (TH2F*)scaleh[1]->Clone();
    ratioscale[0]->Divide(scaleh[1],scaleh[0]);

    TFile* scalesfile2 = TFile::Open("newweights/ZZ_mufhalf.root");
    scaleh[2] = (TH2F*)scalesfile2->Get("Pt_bkg");
    scaleh[2]->Rebin2D(4,4);   scaleh[2]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[2]->Scale(1./scaleh[2]->Integral());
    ratioscale[1] = (TH2F*)scaleh[2]->Clone();
    ratioscale[1]->Divide(scaleh[2],scaleh[0]);

    TFile* scalesfile3 = TFile::Open("newweights/ZZ_murdouble.root");
    scaleh[3] = (TH2F*)scalesfile3->Get("Pt_bkg");
    scaleh[3]->Rebin2D(4,4);   scaleh[3]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[3]->Scale(1./scaleh[3]->Integral());
    ratioscale[2] = (TH2F*)scaleh[3]->Clone();
    ratioscale[2]->Divide(scaleh[3],scaleh[0]);

    TFile* scalesfile4 = TFile::Open("newweights/ZZ_murhalf.root");
    scaleh[4] = (TH2F*)scalesfile4->Get("Pt_bkg");
    scaleh[4]->Rebin2D(4,4);   scaleh[4]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[4]->Scale(1./scaleh[4]->Integral());
    ratioscale[3] = (TH2F*)scaleh[4]->Clone();
    ratioscale[3]->Divide(scaleh[4],scaleh[0]);

    TFile* scalesfile5 = TFile::Open("newweights/ZZ_mufdouble_murdouble.root");
    scaleh[5] = (TH2F*)scalesfile5->Get("Pt_bkg");
    scaleh[5]->Rebin2D(4,4);   scaleh[5]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[5]->Scale(1./scaleh[5]->Integral());
    ratioscale[4] = (TH2F*)scaleh[5]->Clone();
    ratioscale[4]->Divide(scaleh[5],scaleh[0]);

    TFile* scalesfile6 = TFile::Open("newweights/ZZ_mufhalf_murhalf.root");
    scaleh[6] = (TH2F*)scalesfile6->Get("Pt_bkg");
    scaleh[6]->Rebin2D(4,4);   scaleh[6]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[6]->Scale(1./scaleh[6]->Integral());
    ratioscale[5] = (TH2F*)scaleh[6]->Clone();
    ratioscale[5]->Divide(scaleh[6],scaleh[0]);

    TFile* scalesfile7 = TFile::Open("newweights/ZZ_mufhalf_murhalf.root");
    scaleh[7] = (TH2F*)scalesfile7->Get("Pt_bkg");
    scaleh[7]->Rebin2D(4,4);   scaleh[7]->Sumw2();
    // scaleh[0]->GetXaxis()->SetRange(1,80);
    scaleh[7]->Scale(1./scaleh[7]->Integral());
    ratioscale[6] = (TH2F*)scaleh[7]->Clone();
    ratioscale[6]->Divide(scaleh[7],scaleh[0]);  
     
    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	int theBin = scaleh[0]->FindBin(mzz,ptzz);
 	pth->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
	pth1->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*ratioscale[0]->GetBinContent(theBin));
  	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
        // cout << "thisDiff scale-ZZ = " << thisDiff << endl;  
	ptvbf->Fill(mzz,ptVar(ptzz,mzz,overM),wzz*thisDiff); 
      }
    } 
    
    sprintf(nameSyst,"scale-ZZ");
    
  }

  if (whichtype == 0 || whichtype == 3) return;

  // Save and draw histos
  
  char whichSample[4] = "gg";
  if (whichtype == -7) sprintf(whichSample,"zh");
  if (whichtype == -6) sprintf(whichSample,"wh");
  if (whichtype == -4 || whichtype == -5) sprintf(whichSample,"vbf");
  if (whichtype > 0) sprintf(whichSample,"zz");
  
  int uplimit[3];
  int downlimit[3];
  downlimit[0] = pth->GetXaxis()->FindBin(106.);
  uplimit[0]   = pth->GetXaxis()->FindBin(141.);
  downlimit[1] = pth->GetXaxis()->FindBin(213.);
  uplimit[1]   = pth->GetXaxis()->FindBin(288.);
  downlimit[2] = pth->GetXaxis()->FindBin(526.);
  uplimit[2]   = pth->GetXaxis()->FindBin(1257.);

  TH1F* projmass[3];
  projmass[0] = (TH1F*)pth->ProjectionY("projmass0",downlimit[0],uplimit[0]);    // around 126 GeV
  projmass[1] = (TH1F*)pth->ProjectionY("projmass1",downlimit[1],uplimit[1]);   // around 250 GeV
  projmass[2] = (TH1F*)pth->ProjectionY("projmass2",downlimit[2],uplimit[2]);   // around 800 GeV

  TH1F* projmass1[3];
  projmass1[0] = (TH1F*)pth1->ProjectionY("projmass10",downlimit[0],uplimit[0]);    // around 126 GeV
  projmass1[1] = (TH1F*)pth1->ProjectionY("projmass11",downlimit[1],uplimit[1]);   // around 250 GeV
  projmass1[2] = (TH1F*)pth1->ProjectionY("projmass12",downlimit[2],uplimit[2]);   // around 800 GeV

  TH1F* projmass2[3];
  projmass2[0] = (TH1F*)pth2->ProjectionY("projmass20",downlimit[0],uplimit[0]);    // around 126 GeV
  projmass2[1] = (TH1F*)pth2->ProjectionY("projmass21",downlimit[1],uplimit[1]);   // around 250 GeV
  projmass2[2] = (TH1F*)pth2->ProjectionY("projmass22",downlimit[2],uplimit[2]);   // around 800 GeV

  TH1F* projmassvbf[3];
  projmassvbf[0] = (TH1F*)ptvbf->ProjectionY("projmassvbf0",downlimit[0],uplimit[0]);    // around 126 GeV
  projmassvbf[1] = (TH1F*)ptvbf->ProjectionY("projmassvbf1",downlimit[1],uplimit[1]);   // around 250 GeV
  projmassvbf[2] = (TH1F*)ptvbf->ProjectionY("projmassvbf2",downlimit[2],uplimit[2]);   // around 800 GeV
  
  for (int k = 0; k < 3; k++) {
    
    projmass[k]->Scale(1./projmass[k]->Integral());
    projmass1[k]->Scale(1./projmass1[k]->Integral());
    projmass2[k]->Scale(1./projmass2[k]->Integral());
    projmassvbf[k]->Scale(1./projmassvbf[k]->Integral());

    if (overM == 0) sprintf(nameFile,"p_{T} [GeV/c]");
    else if (overM == 1) sprintf(nameFile,"p_{T}/m_{4l}");
    else sprintf(nameFile,"ln(p_{T}/m_{4l})");

    projmass[k]->SetMarkerColor(1);    projmass[k]->SetLineColor(1);    projmass[k]->SetLineWidth(2);
    projmass1[k]->SetMarkerColor(2);   projmass1[k]->SetLineColor(2);   projmass1[k]->SetLineWidth(2);
    projmass2[k]->SetMarkerColor(4);   projmass2[k]->SetLineColor(4);   projmass2[k]->SetLineWidth(2);
    projmassvbf[k]->SetMarkerColor(6);   projmassvbf[k]->SetLineColor(6);   projmassvbf[k]->SetLineWidth(2);
  			    			   
    can.cd(1);
    gPad->SetBottomMargin(0.0);
    if (whichtype == -3) gPad->SetLogy();
    projmass2[k]->GetXaxis()->SetLabelColor(kWhite);
    projmass2[k]->GetYaxis()->SetTitle("Normalized entries");
    projmass2[k]->Draw();
    if (whichtype != -3 && whichtype != -6 && whichtype != -7 && whichtype != 4) projmass1[k]->Draw("SAME");
    projmass[k]->Draw("SAME");
    // if (whichtype == -5 || whichtype == -4 || whichtype == 6 || whichtype == 7) projmassvbf[k]->Draw("SAME");
  
    TH1F* diffpt = (TH1F*)projmass[k]->Clone();
    TH1F* pullpt = (TH1F*)projmass[k]->Clone();
    diffpt->Add(projmassvbf[k],projmass[k],1,-1);
    pullpt->Divide(diffpt,projmass[k]);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-1.);
    pullpt->SetMaximum(1.);
    pullpt->GetXaxis()->SetTitle(nameFile);
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
  
    char mass[10];
    if (k == 0) sprintf(mass,"126");
    if (k == 1) sprintf(mass,"250");
    if (k == 2) sprintf(mass,"800");

    sprintf(nameFile,"newfigs/%s_syst%s_%sGeV.gif",LcasePt,nameSyst,mass);
    if (!also7TeV || (whichtype == 5 && also7TeV) ) can.SaveAs(nameFile);
    sprintf(nameFile,"newfigs/%s_syst%s_%sGeV.pdf",LcasePt,nameSyst,mass);
    if (!also7TeV || (whichtype == 5 && also7TeV) ) can.SaveAs(nameFile);
    
  }

  sprintf(nameFile,"selRootFiles/%s_%s_TEMPL_8TeV.root",UcasePt,whichSample);
  if (also7TeV) sprintf(nameFile,"selRootFiles/%s_%s_TEMPL_7TeV.root",UcasePt,whichSample);
  TFile f1(nameFile,"UPDATE");
  can2.cd();   gPad->SetLogx();   
  sprintf(nameFile,"%sH_%s",LcasePt,nameSyst);
     
  adjustHistogram(pth);
  adjustHistogram(pth1);
  adjustHistogram(pth2);
  adjustHistogram(ptvbf); 

  // cout << whichtype << " " << nameFile << endl;
  if (whichtype == -4 || whichtype == -5 || whichtype == 6 || whichtype == 7) {
    ptvbf->SetName(nameFile);
    ptvbf->Write();
    ptvbf->Draw("COLZ");
  } else if (whichtype == -3 || whichtype == 3 || whichtype == 4 || whichtype == 5) {
    pth->SetName(nameFile);
    pth->Write();
    pth->Draw("COLZ");
  } else if (whichtype == -2 || whichtype == -6 || whichtype == -7) {
    pth2->SetName(nameFile);
    pth2->Write();
    pth2->Draw("COLZ");
  }
  /* sprintf(nameFile,"%sH_Default",LcasePt);
  TH1F* def = (TH1F*)f1.Get(nameFile);
  if (whichtype == -4 || whichtype == -5 || whichtype == 6 || whichtype == 7) evalBinMigration(def,ptvbf,nameSyst,true);
  else if (whichtype == -3 || whichtype == 3 ||  whichtype == 4 || whichtype == 5) evalBinMigration(def,pth,nameSyst,true);
  else if (whichtype == -2 || whichtype == -6 || whichtype == -7) evalBinMigration(def,pth2,nameSyst,true); */
  f1.Close();
  sprintf(nameFile,"newfigs/%s_%s%sTemplate_%dTeV.gif",LcasePt,whichSample,nameSyst,sqrts);
  can2.SaveAs(nameFile);
  sprintf(nameFile,"newfigs/%s_%s%sTemplate_%dTeV.pdf",LcasePt,whichSample,nameSyst,sqrts);
  can2.SaveAs(nameFile);
  
  return;
}

