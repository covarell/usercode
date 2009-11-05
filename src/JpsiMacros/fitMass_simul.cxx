// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

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
#include "RooSimWSTool.h"

using namespace RooFit;

void defineBackground(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // BKG: first and second order polynomials
  RooRealVar coefPol1("coefPol1","linear coefficient of bkg PDF",-0.05,-150.,150.);
  RooRealVar coefPol2("coefPol2","quadratic coefficient of bkg PDF",0.1,-1.,1.);

  RooRealVar CcoefPol1("CcoefPol1","linear coefficient of bkg PDF",-0.05,-1500.,1500.);
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
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.1);
  RooExponential expFunct("expFunct","expFunct",*JpsiMass,coefExp); 

  ws->import(RooArgSet(PolFunct1,PolFunct2,CPolFunct1,CPolFunct2,expFunct));

  return;
}

void defineSignal(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // SIGNAL: 2-Gaussians
  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.1,3.05,3.15);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.02,0.001,0.2);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.1,3.05,3.15);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.03,0.001,0.2);

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
  RooRealVar enne("enne","n of CB",3.,0.,10.);

  RooCBShape sigCB("sigCB","Signal CB PDF",*JpsiMass,meanSig1,sigmaSig1,alpha,enne);

  RooAddPdf sigCBGauss("sigCBGauss","Signal CB+Gauss PDF",sigCB,signalG2,coeffGauss);

  ws->import(RooArgSet(sigPDF,sigPDFOneMean,sigCBGauss),RecycleConflictNodes());

  return;
}

void setRanges(RooWorkspace *ws){

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;

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

void drawResults(RooWorkspace *ws, const bool isGG){

  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totalPDF = ws->pdf("TOTPDF_sim");
  RooAbsPdf *background = ws->pdf("expFunct");
  RooAbsPdf *signal = ws->pdf("sigCBGauss");
  RooCategory *JpsiPtType = ws->cat("JpsiPtType");
  RooCategory *JpsiEtaType = ws->cat("JpsiEtaType");
  RooRealVar *JpsiMass = ws->var("JpsiMass");

  string reducestr;

  Int_t nbinspt = 11;
  Int_t nbinseta = 2;

  TCanvas c1;
  c1.Divide(3,4);

  TCanvas c2;
  c2.Divide(3,4);

  for(int i = 1;i<=nbinspt;i++){
    for(int j = 1;j<=nbinseta;j++){

      RooPlot *mframe = JpsiMass->frame();

      ostringstream ost;
      ost << i << "th pT bin and in the " << j << "th eta bin";

      if(isGG) reducestr = "Mass fit for glb-glb muons in the " + ost.str();
      else reducestr = "Mass fit for glb-trk muons in the " + ost.str();

      mframe->SetTitle(reducestr.c_str());

      string cutstring;

      ostringstream ostdrei;
      ostdrei << "P" << i;
      ostringstream ostvier;
      ostvier << "E" << j;

      string catstring_pt = ostdrei.str();
      string catstring_eta = ostvier.str();

      if(isGG) cutstring = "JpsiType == JpsiType::GG && JpsiPtType == JpsiPtType::" + catstring_pt + " && JpsiEtaType == JpsiEtaType::" + catstring_eta;
      else cutstring = "JpsiType == JpsiType::GT && JpsiPtType == JpsiPtType::" + catstring_pt + " && JpsiEtaType == JpsiEtaType::" + catstring_eta;

      data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut(cutstring.c_str()));

      JpsiPtType->setLabel(catstring_pt.c_str());
      JpsiEtaType->setLabel(catstring_eta.c_str());

      string bkgname = "expFunct_{" + catstring_pt + ";" + catstring_eta + "}";

      string totPDFname = "totPDF_{" + catstring_pt + ";" + catstring_eta + "}";

      ws->pdf(totPDFname.c_str())->plotOn(mframe,Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));
      ws->pdf(totPDFname.c_str())->plotOn(mframe,DrawOption("F"),FillColor(kGreen),Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));
      ws->pdf(totPDFname.c_str())->plotOn(mframe,Components(*(ws->pdf(bkgname.c_str()))),DrawOption("F"),FillColor(kRed),Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));

      data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut(cutstring.c_str()));

      if(j == 1) {c1.cd(i);mframe->Draw();}
      if(j == 2) {c2.cd(i);mframe->Draw();}

    }
  }

  if(isGG) reducestr = "GGmassfit_sim.gif";
  else reducestr = "GTmassfit_sim.gif";
  c1.SaveAs(reducestr.c_str());

  if(isGG) reducestr = "GGmassfit_sim2.gif";
  else reducestr = "GTmassfit_sim2.gif";
  c2.SaveAs(reducestr.c_str());

  return;
}

void printResults(RooWorkspace *ws){

  Int_t nbinspt = 11;
  Int_t nbinseta = 2;

  cout << endl << "GG J/psi yields:" << endl;

  for(int i = 1;i<=nbinspt;i++){
    for(int j = 1;j<=nbinseta;j++){

      string cutstring;

      ostringstream ostdrei;
      ostdrei << "P" << i;
      ostringstream ostvier;
      ostvier << "E" << j;

      string catstring_pt = ostdrei.str();
      string catstring_eta = ostvier.str();

      string Nsigname = "NSig_{" + catstring_pt + ";" + catstring_eta + "}";
      string coeffGaussname = "coeffGauss_{" + catstring_pt + ";" + catstring_eta + "}";

      const double Nsig   = ws->var(Nsigname.c_str())->getVal();
      const double errSig = ws->var(Nsigname.c_str())->getError();
      const double coeffGauss = ws->var(coeffGaussname.c_str())->getVal();
      const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
      const double sigmaSig2 = ws->var("sigmaSig2")->getVal();

      const double resol = (coeffGauss*coeffGauss*sigmaSig1 + (1-coeffGauss)*(1-coeffGauss)*sigmaSig2) / (coeffGauss*coeffGauss + (1-coeffGauss)*(1-coeffGauss)); 

      cout << "Fit for bin (" << i << "," << j << "): " << Nsig << " +/- " << errSig << endl;
      cout << "Resolution : " << resol*1000. << " MeV" << endl;
    }
  }

  return;
}

int main(int argc, char* argv[]) {

  /// LIMITS ///

  gROOT->SetStyle("Plain");

  char *filename;

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
      
      }
    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *reddata1 = (RooDataSet*)fIn.Get("data");

  reddata1->setWeightVar("MCweight");

  ws->import(*reddata1);

  setRanges(ws);

  ws->var("JpsiMass")->setBins(30);

  RooDataHist *reddata = new RooDataHist("reddata","reddata",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType")),*(ws->cat("JpsiType"))),*reddata1);

  //DEFINE SIGNAL AND BACKGROUND
  defineBackground(ws);
  defineSignal(ws);

  RooRealVar NSig("NSig","Number of signal events",5000.,1.,10000000.);
  RooRealVar NBkg("NBkg","Number of background events",2000.,1.,10000000.);

  // Total PDF (signal CB+Gauss)
  RooAddPdf totPDF("totPDF","Total pdf",RooArgList(*(ws->pdf("sigCBGauss")),*(ws->pdf("expFunct"))),RooArgList(NSig,NBkg));

  ws->import(totPDF);

  //Create the simultaneous PDF for the pt bins fit
  RooSimWSTool wst(*ws);
  wst.build("TOTPDF_sim","totPDF",SplitParam("coefExp,coeffGauss,NSig,NBkg","JpsiPtType,JpsiEtaType")) ;

  //RooSimWSTool::MultiBuildConfig mbc("JpsiEtaType");
  //mbc.addPdf("E1","totPDF",SplitParam("coefExp,coeffGauss,,NSig,NBkg","JpsiPtType"));
  //mbc.addPdf("E2","totPDF",SplitParam("coefExp,coeffGauss,,NSig,NBkg","JpsiPtType"));
  //mbc.restrictBuild("JpsiPtType","P3,P4,P5,P6,P7,P8,P9,P10,P11");
  //wst.build("TOTPDF_sim",mbc);

  RooDataHist *GGdata = (RooDataHist*)reddata->reduce("JpsiType == JpsiType::GG");
  ws->pdf("TOTPDF_sim")->fitTo(*GGdata,Extended(1),Save(1),Minos(0));

  drawResults(ws,true);

  printResults(ws);

  // fix some parameters 
  ws->var("alpha")->setConstant(kTRUE); 
  ws->var("enne")->setConstant(kTRUE); 

  //GT

  RooDataHist *GTdata = (RooDataHist*)reddata->reduce("JpsiType == JpsiType::GT");
  ws->pdf("TOTPDF_sim")->fitTo(*GTdata,Extended(1),Save(1),Minos(0));

  drawResults(ws,false);

  printResults(ws);

  return 1;
}
