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
#include "PDFs/RooModifTsallis.h"
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

Double_t zzFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 182.) fitval = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  else fitval = par[3]*x[0] + par[4]*x[0]*x[0] - par[3]*182. - par[4]*33124. + par[0] + par[1]*182. + par[2]*33124.;
  return fitval;
}

/* Double_t vbfFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 300.) fitval = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  else fitval = par[3]*x[0] + par[4]*x[0]*x[0] - par[3]*300. - par[4]*90000. + par[0] + par[1]*300. + par[2]*90000.;
  return fitval;
  }*/

Double_t fexpFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 250.) fitval = 0.005;
  else fitval = par[0]*x[0] + par[1]*x[0]*x[0] + par[2]*x[0]*x[0]*x[0] - par[0]*250. - par[1]*62500. - par[2]*15625000. + 0.005;
  return fitval;
}

Double_t BBFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 400.) fitval = par[0];
  else if (x[0] >= 700.) fitval = par[1];
  else fitval = par[0] - (par[1]-par[0])*400./300. + (par[1]-par[0])*x[0]/300.;
  return fitval;
}

Double_t vhFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 140.) fitval = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
  else fitval = par[0] + par[1]*140. + par[2]*19600. + par[3]*2744000.;
  return fitval;
}

//=======================================================================

TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,int LHCsqrts= 7,bool isLowMass=true, TString systName = "Default", bool down = false){  
 
  char fileName[200];

  static const int ptmbins = 20;
  double ptmbinLimits[ptmbins+1] = {0.,.1,.2,.3,.4,.5,.6,.7,
                                    .8,.9,1.0,1.1,1.2,1.3,1.4,1.5,
                                    1.75,2.0,2.25,2.5,4.};

  TH2F* bkgHist;
  if(!isLowMass)
    // Fine binning (constant)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,500,0,4);
    // Coarse binning (constant)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,50,0,3);
    // Coarse binning (variable)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,ptmbins,ptmbinLimits);
  else
    // Fine binning (constant)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,500,0,4);
    // Coarse binning (constant)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,50,0,3);
    // Coarse binning (variable)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,ptmbins,ptmbinLimits);

  bkgHist->Sumw2();

  if (sampleIndex < 0 || sampleIndex > nsamp) {
    cout << "Samples are numbered from 0 to " << nsamp-1 << endl;
    return bkgHist;
  }

  // functions of mZZ
  TF1* mfunc;
  TF1* nfunc;
  TF1* n2func;
  TF1* bbfunc;
  TF1* bb2func;
  TF1* Tfunc;
  TF1* fexpfunc;

  if (sampleIndex == 0) {  // gg
    mfunc = new TF1("mfunc",zzFunc,100.,1600.,5);
    nfunc = new TF1("nfunc","pol0");
    n2func = new TF1("n2func","pol0");
    bbfunc = new TF1("bbfunc",zzFunc,100.,1600.,5);
    bb2func = new TF1("bb2func",zzFunc,100.,1600.,5);
    Tfunc = new TF1("Tfunc",zzFunc,100.,1600.,5);
    fexpfunc = new TF1("fexpfunc","pol0");
  } else if (sampleIndex == 1) {  // vbf
    mfunc = new TF1("mfunc","pol3");
    nfunc = new TF1("nfunc","pol0");
    n2func = new TF1("n2func","pol0");
    bbfunc = new TF1("bbfunc",BBFunc,100.,1600.,2);
    bb2func = new TF1("bb2func",BBFunc,100.,1600.,2);
    Tfunc = new TF1("Tfunc","pol4");
    fexpfunc = new TF1("fexpfunc",fexpFunc,100.,1600.,3);
  } else if (sampleIndex == 4) {  // zz
    mfunc = new TF1("mfunc",zzFunc,100.,1600.,5);
    nfunc = new TF1("nfunc","pol0");
    n2func = new TF1("n2func","pol0");
    bbfunc = new TF1("bbfunc",zzFunc,100.,1600.,5);
    bb2func = new TF1("bb2func","pol0");
    Tfunc = new TF1("Tfunc","pol0");
    fexpfunc = new TF1("fexpfunc","pol0");
  } else {  // wh/zh
    mfunc = new TF1("mfunc",vhFunc,100.,1600.,4);
    nfunc = new TF1("nfunc","pol0");
    n2func = new TF1("n2func","pol0");
    bbfunc = new TF1("bbfunc","pol0");
    bb2func = new TF1("bb2func","pol0");
    Tfunc = new TF1("Tfunc",vhFunc,100.,1600.,4);
    fexpfunc = new TF1("fexpfunc","pol0");
  }
    
  
  //read parameters
  char thePar[10];
  char thePar2[3];
  char equalS[1];
  float fitted;
  int ii;

  if (systName == "OneSyst") {
    if (down) sprintf(fileName,"tsallisPars/allParamsDown_oneSyst_%s_%dTeV.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
    else sprintf(fileName,"tsallisPars/allParamsUp_oneSyst_%s_%dTeV.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
  } else sprintf(fileName,"tsallisPars/allParams_%s_%dTeV.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
  cout << fileName << endl;
  ifstream theFile(fileName);
  
  while (theFile >> thePar >> thePar2 >> equalS >> fitted) {
    sscanf(thePar2,"[%d]",&ii);
    if (sampleIndex == 3) cout << thePar << " [" << ii << "] = " << fitted << " " << endl;
    if (!strcmp(thePar,"bb") ) bbfunc->SetParameter(ii,fitted);
    if (!strcmp(thePar,"bbdue") ) bb2func->SetParameter(ii,fitted);
    if (!strcmp(thePar,"n") ) nfunc->SetParameter(ii,fitted);
    if (!strcmp(thePar,"ndue") ) n2func->SetParameter(ii,fitted);
    if (!strcmp(thePar,"m") ) mfunc->SetParameter(ii,fitted);
    if (!strcmp(thePar,"T") ) Tfunc->SetParameter(ii,fitted);
    if (!strcmp(thePar,"fexp") ) fexpfunc->SetParameter(ii,fitted);
  }

  // fill histogram
  RooRealVar* ptoverm = new RooRealVar("ptoverm","p_{T}/M^{4l}",0.,4.,"GeV/c");

  RooRealVar mup("mup","emme", 1.,-1000000.,1000000.);
  RooRealVar nup("nup","enne", 0.93, -1000000.,1000000.);
  RooRealVar n2up("n2up","enne2", 0.75, -1000000.,1000000.);
  RooRealVar bbup("bbup","bibi",0.02, -1000000.,1000000.);
  RooRealVar Tup("Tup","tti",0.2,-1000000.,1000000.);
  RooRealVar bb2up("bb2up","bibi2",0.02, -1000000.,1000000.);
  RooRealVar fexpup("fexpup","f_exp",0.02, -1000000.,1000000.);
 
  RooModifTsallis* rtup = new RooModifTsallis("rtup","rtup",*ptoverm,
					    mup,nup,n2up,bbup,bb2up,Tup,fexpup);

  // default
  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();
    
  // if (systName == "Default") {
    
    for (Int_t j=1; j<=nXbins; j++) {      
      float mzz = bkgHist->GetXaxis()->GetBinCenter(j);
      mup.setVal(mfunc->Eval(mzz));
      nup.setVal(nfunc->Eval(mzz));
      n2up.setVal(n2func->Eval(mzz));
      bbup.setVal(bbfunc->Eval(mzz));
      bb2up.setVal(bb2func->Eval(mzz));
      fexpup.setVal(fexpfunc->Eval(mzz));
      Tup.setVal(Tfunc->Eval(mzz));
      // cout << mzz << " " << mfunc->Eval(mzz) << " " << nfunc->Eval(mzz) << endl;
      for (Int_t i=1; i<=nYbins; i++) {
	ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
	bkgHist->SetBinContent(j,i,rtup->getVal());
      }
    } 

  
  // normalize slices

  double norm;
  TH1F* tempProj;
  
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
  
  return bkgHist;
  
}

//=======================================================================

void makeTemplate(TString channel="4mu",bool isLowMass=true){

  char fileName[200];

  for (Int_t lhc=8; lhc<9; lhc++) {
  // for (Int_t lhc=7; lhc<9; lhc++) {
    for (Int_t k=0; k<nsamp; k++) {
     
      TString lhcs = "7";
      if (lhc==8) lhcs="8";

      TFile* ftemp = new TFile(destDir + "PtoverMComplicated_mZZ_" + dataFileNames[k] + "_" + channel + "_"+ lhcs + "TeV.root","RECREATE");

      cout << "*** Now filling: " << dataFileNames[k] << ", Default template" << endl;
      TH2F* h_mzzD = fillTemplate(channel,k,lhc,isLowMass);
     
      ftemp->cd();
  
      h_mzzD->Write("h_Ptmzz_mzz");
      cout << " " << h_mzzD->GetBinContent(100,10) << endl;
      
      cout << "*** Now filling: " << dataFileNames[k] << ", OneSyst templates" << endl;
      
      TH2F* h_mzzDup = fillTemplate(channel,k,lhc,isLowMass,"OneSyst",false);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_up");
      h_mzzDup->Write(fileName);
      cout << " " << h_mzzDup->GetBinContent(100,10) << endl;
      
      TH2F* h_mzzDdown = fillTemplate(channel,k,lhc,isLowMass,"OneSyst",true);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_down");
      h_mzzDdown->Write(fileName);
      cout << " " << h_mzzDdown->GetBinContent(100,10) << endl;

      ftemp->Close();

    }
  }

}

//=======================================================================

void generateComplicatedPtTemplates(){

  bool isLowMass = false;

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

  makeTemplate("4mu",isLowMass);
  makeTemplate("4e",isLowMass);
  makeTemplate("2e2mu",isLowMass);

}



