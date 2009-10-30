// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

using namespace RooFit;

void defineBackground(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // BKG: first and second order polynomials
  RooRealVar coefPol1("coefPol1","linear coefficient of bkg PDF",-0.05,-15000.,15000.);
  RooRealVar coefPol2("coefPol2","quadratic coefficient of bkg PDF",0.1,-1.,1.);

  RooRealVar CcoefPol1("CcoefPol1","linear coefficient of bkg PDF",-0.05,-15000.,15000.);
  RooRealVar CcoefPol2("CcoefPol2","quadratic coefficient of bkg PDF",0.1,-1.,1.);

  RooPolynomial PolFunct1("PolFunct","PolFunct1",*JpsiMass,coefPol1);
  RooPolynomial PolFunct2("PolFunct2","PolFunct2",*JpsiMass,RooArgList(coefPol1,coefPol2));

  RooPolynomial CPolFunct1("CPolFunct","CPolFunct1",*JpsiMass,CcoefPol1);
  RooPolynomial CPolFunct2("CPolFunct2","CPolFunct2",*JpsiMass,RooArgList(CcoefPol1,CcoefPol2));
  
  coefPol2.setVal(0.0); 
  coefPol2.setConstant(kTRUE); 
  CcoefPol2.setVal(0.0); 
  CcoefPol2.setConstant(kTRUE); 

  // BKG : exponential
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.);
  RooExponential expFunct("expFunct","expFunct",*JpsiMass,coefExp); 

  ws->import(RooArgSet(PolFunct1,PolFunct2,CPolFunct1,CPolFunct2,expFunct));

  return;
}

void defineSignal(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // SIGNAL: 2-Gaussians
  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.1,3.05,3.15);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.02,0.,0.2);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.1,3.05,3.15);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.03,0.,0.2);

  // Different mean
  RooGaussian signalG1("signalG1","Signal PDF 1",*JpsiMass,meanSig1,sigmaSig1);
  RooGaussian signalG2("signalG2","Signal PDF 2",*JpsiMass,meanSig2,sigmaSig2);
  // Same mean
  RooGaussian signalG2OneMean("signalG2OneMean","Signal PDF 2",*JpsiMass,meanSig1,sigmaSig2);

  RooRealVar coeffGauss("coeffGauss","Relative norm of the two signal gaussians",0.42,0.,1.);

  // Different mean
  RooAddPdf sigPDF("sigPDF","Total signal pdf",signalG1,signalG2,coeffGauss);
  // Same mean
  RooAddPdf sigPDFOneMean("sigPDFOneMean","Total signal pdf",signalG1,signalG2OneMean,coeffGauss);

  // SIGNAL: Crystal Ball shape
  RooRealVar alpha("alpha","#alpha of CB",0.96,0.,5.);
  RooRealVar enne("enne","n of CB",3.,0.,8.);

  RooCBShape sigCB("sigCB","Signal CB PDF",*JpsiMass,meanSig1,sigmaSig1,alpha,enne);

  RooAddPdf sigCBGauss("sigCBGauss","Signal CB+Gauss PDF",sigCB,signalG2,coeffGauss);

  ws->import(RooArgSet(sigPDF,sigPDFOneMean,sigCBGauss),RecycleConflictNodes());

  return;
}

void getrange(string &varRange, float *varmin, float *varmax){
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
    cout << "pT range not valid!" << endl;
    assert(0);
  }

 return;
}

void prefitSideband(RooWorkspace *ws, bool isGG){

  RooDataSet *tmpdata;
  if(isGG) tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GG");
  else tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GT");

  if(isGG) ws->pdf("PolFunct2")->fitTo(*tmpdata,Range("left,right"));
  else ws->pdf("PolFunct2")->fitTo(*tmpdata,Range("left,right"));

  ws->var("CcoefPol1")->setVal(ws->var("coefPol1")->getVal());
  ws->var("CcoefPol1")->setConstant(kTRUE);
  ws->var("CcoefPol2")->setVal(ws->var("coefPol2")->getVal());
  ws->var("CcoefPol2")->setConstant(kTRUE);

  // float theCurrentVal = coefPol1.getVal();
  // coefPol1.setRange(theCurrentVal-0.01,theCurrentVal+0.01);

  return;
}

void setRanges(RooWorkspace *ws){

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;
  //const float JpsiCtMin = -1.0;
  //const float JpsiCtMax = 5.0;

  ws->var("JpsiMass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("JpsiMass")->setRange("left",JpsiMassMin,2.9);
  ws->var("JpsiMass")->setRange("right",3.3,JpsiMassMax);

  ws->cat("JpsiType")->setRange("glbglb","GG");
  ws->cat("JpsiType")->setRange("glbtrk","GT");

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

void drawResults(RooWorkspace *ws, const bool isGG, const string prange, const string etarange){

  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDFCBGauss = ws->pdf("totPDFCBGauss");
  RooAbsPdf *CPolFunct2 = ws->pdf("CPolFunct2");
  RooAbsPdf *sigCBGauss = ws->pdf("sigCBGauss");

  char reducestr[200];

  RooPlot *mframe = ws->var("JpsiMass")->frame();

  if(isGG) sprintf(reducestr,"Mass fit for glb-glb muons p_{T} = %s GeV,   #eta = %s",prange.c_str(),etarange.c_str());
  else sprintf(reducestr,"Mass fit for glb-trk muons p_{T} = %s GeV,   #eta = %s",prange.c_str(),etarange.c_str());

  mframe->SetTitle(reducestr);

  if(isGG) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  // PolFunct2.plotOn(GGmframe,Range("left,right"));
  totPDFCBGauss->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDFCBGauss->plotOn(mframe,Components(RooArgSet(*CPolFunct2,*sigCBGauss)),DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDFCBGauss->plotOn(mframe,Components(*CPolFunct2),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
  // totPDFCBGaussExp.plotOn(GGmframe);
  // totPDFCBGaussExp.plotOn(GGmframe,Components(RooArgSet(expFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  // totPDFCBGaussExp.plotOn(GGmframe,Components(expFunct),DrawOption("F"),FillColor(kRed));

  if(isGG) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  TCanvas c1;
  c1.cd();mframe->Draw();
  if(isGG) sprintf(reducestr,"GGmassfit_pT%s_eta%s.gif",prange.c_str(),etarange.c_str());
  else sprintf(reducestr,"GTmassfit_pT%s_eta%s.gif",prange.c_str(),etarange.c_str());
  c1.SaveAs(reducestr);

  return;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol){

  Nsig   = ws->var("NSig")->getVal();
  errSig = ws->var("NSig")->getError();
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  const double sigmaSig2 = ws->var("sigmaSig2")->getVal();

  resol = (coeffGauss*coeffGauss*sigmaSig1 + (1-coeffGauss)*(1-coeffGauss)*sigmaSig2) / (coeffGauss*coeffGauss + (1-coeffGauss)*(1-coeffGauss)); 

  return;
}

int main(int argc, char* argv[]) {

  /// LIMITS ///

  gROOT->SetStyle("Plain");

  char *filename;
  string prange;
  string etarange;
  bool sidebandPrefit = false;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];
    
    switch(pchar[0]){
      
    case '-':{
      
      switch(pchar[1]){
	
      case 'f':{
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
      
      case 'p':{
	prange = argv[i+1];
	cout << "Range for pT is " << prange << " GeV/c" << endl;
        break;
      }
       
      case 'e':{
        etarange = argv[i+1];
        cout << "Range for |eta| is " << etarange << endl;
        break;
      }

      case 's':{
        sidebandPrefit = true;
        cout << "Sideband pre-fitting activated" << endl;
        break;
      }
      }
    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");
  
  float pmin, pmax; 
  float etamin, etamax;

  getrange(prange,&pmin,&pmax);
  getrange(etarange,&etamin,&etamax);

  char reducestr[200];
  sprintf(reducestr,"JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f", pmax,pmin,etamax,etamin);
  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);

  ws->import(*reddata);

  setRanges(ws);

  //DEFINE SIGNAL AND BACKGROUND
  defineBackground(ws);
  defineSignal(ws);

  RooRealVar NSig("NSig","Number of signal events",5000.,10.,10000000.);
  RooRealVar NBkg("NBkg","Number of background events",2000.,10.,10000000.);

  // Total PDF (signal CB+Gauss)
  RooAddPdf totPDFCBGauss("totPDFCBGauss","Total pdf",RooArgList(*(ws->pdf("sigCBGauss")),*(ws->pdf("CPolFunct2"))),RooArgList(NSig,NBkg));

  ws->import(totPDFCBGauss);

  /*
  // Total PDF 
  RooAddPdf totPDF("totPDF","Total pdf",RooArgList(sigPDF,CPolFunct2),RooArgList(NSig,NBkg));
  // Total PDF (exponential background)
  RooAddPdf totPDFExp("totPDFexp","Total pdf",RooArgList(sigPDF,expFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean)
  RooAddPdf totPDFOneMean("totPDFOneMean","Total pdf",RooArgList(sigPDFOneMean,CPolFunct2),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean and exponential background)
  RooAddPdf totPDFExpOneMean("totPDFexpOneMean","Total pdf",RooArgList(sigPDFOneMean,expFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean and quadratic background)
  RooAddPdf totPDF2Pol("totPDF2Pol","Total pdf",RooArgList(sigPDFOneMean,CPolFunct2),RooArgList(NSig,NBkg));
  // Total PDF (signal CB)
  RooAddPdf totPDFCB("totPDFCB","Total pdf",RooArgList(sigCB,CPolFunct2),RooArgList(NSig,NBkg));
  // Total PDF (signal CB and exponential background)
  RooAddPdf totPDFExpCB("totPDFExpCB","Total pdf",RooArgList(sigCB,expFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal CB+Gauss and exponential background)
  RooAddPdf totPDFExpCBGauss("totPDFExpCBGauss","Total pdf",RooArgList(sigCBGauss,expFunct),RooArgList(NSig,NBkg));
  */

  //GG
  // fix some parameters 
  // alpha.setConstant(kTRUE); 
  // enne.setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,true);

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  ws->pdf("totPDFCBGauss")->fitTo(*GGdata,Extended(1),Save(1),Minos(0));

  drawResults(ws,true,prange,etarange);

  double NSigGG, errSigGG,resolGG;
  printResults(ws,NSigGG,errSigGG,resolGG);

  // fix some parameters 
  ws->var("alpha")->setConstant(kTRUE); 
  ws->var("enne")->setConstant(kTRUE); 
  // CcoefPol1.setVal(-1.0);
  // CcoefPol2.setConstant(kFALSE);

  //GT
  if (sidebandPrefit) prefitSideband(ws,false);

  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  ws->pdf("totPDFCBGauss")->fitTo(*GTdata,Extended(1),Save(1),Minos(0));
  // totPDFExpCBGauss.fitTo(*GTdata,Extended(1),Save(1),Minos(0));

  drawResults(ws,false,prange,etarange);

  double NSigGT, errSigGT,resolGT;
  printResults(ws,NSigGT,errSigGT,resolGT);

  RooDataSet *GGdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG && (MCType == MCType::PR || MCType == MCType::NP)");

  RooDataSet *GTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT && (MCType == MCType::PR || MCType == MCType::NP)");

  cout << endl << "pT = " << prange << " GeV; |eta| = " << etarange << endl;
  cout << endl << "GG J/psi yields:                 GT J/psi yields:" << endl;
  cout << "True MC : " << GGdataTr->numEntries(true) << "                   True MC : " << GTdataTr->numEntries(true) << endl;
  cout << "Fit : " << NSigGG << " +/- " << errSigGG << "        Fit : " << NSigGT << " +/- " << errSigGT << endl;
  cout << "Resolution : " << resolGG*1000. << " MeV         Resolution : " << resolGT*1000. << " MeV" << endl; 

  return 1;
}
