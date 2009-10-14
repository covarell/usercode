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

using namespace RooFit;

int main(int argc, char* argv[]) {

  /// LIMITS ///

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;
  const float JpsiCtMin = -1.0;
  const float JpsiCtMax = 5.0;

  gROOT->SetStyle("Plain");

  char *filename;
  char *prange;
  char *etarange;
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

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");
  
  float pmin, pmax; 
  float etamin, etamax;
  
  if (sscanf(prange, "%f-%f", &pmin, &pmax) == 0) {
    cout << "pT range not valid!" << endl;
    return 0;
  }

  if (sscanf(etarange, "%f-%f", &etamin, &etamax) == 0) {
    cout << "eta range not valid!" << endl;
    return 0;
  }

  char reducestr[200];
  sprintf(reducestr,"JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f", pmax,pmin,etamax,etamin);
  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);

  RooRealVar JpsiMass("JpsiMass","#mu^{+}#mu^{-} mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  JpsiMass.setRange("all",JpsiMassMin,JpsiMassMax);
  JpsiMass.setRange("left",JpsiMassMin,2.9);
  JpsiMass.setRange("right",3.3,JpsiMassMax);

  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi pt",0.,200.,"GeV/c");
  RooRealVar* JpsiEta = new RooRealVar("JpsiEta","J/psi eta",-2.7,2.7);
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,25.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  // Fit functions and parameters

  // BKG: first and second order polynomials
  RooRealVar coefPol1("coefPol1","linear coefficient of bkg PDF",-0.5,-15000.,15000.);
  RooRealVar coefPol2("coefPol2","quadratic coefficient of bkg PDF",0.,-1.,1.);
  RooRealVar CcoefPol1("CcoefPol1","linear coefficient of bkg PDF",-0.5,-15000.,15000.);
  RooRealVar CcoefPol2("CcoefPol2","quadratic coefficient of bkg PDF",0.,-1.,1.);
  RooPolynomial PolFunct1("PolFunct","PolFunct1",JpsiMass,coefPol1);
  RooPolynomial PolFunct2("PolFunct2","PolFunct2",JpsiMass,RooArgList(coefPol1,coefPol2));
  RooPolynomial CPolFunct1("CPolFunct","CPolFunct1",JpsiMass,CcoefPol1);
  RooPolynomial CPolFunct2("CPolFunct2","CPolFunct2",JpsiMass,RooArgList(CcoefPol1,CcoefPol2));
  
  coefPol2.setVal(0.0); 
  coefPol2.setConstant(kTRUE); 
  CcoefPol2.setVal(0.0); 
  CcoefPol2.setConstant(kTRUE); 

  // BKG : exponential
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.);
  RooExponential expFunct("expFunct","expFunct",JpsiMass,coefExp); 

  // SIGNAL: 2-Gaussians
  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.1,3.05,3.15);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.02,0.,0.2);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.1,3.05,3.15);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.03,0.,0.2);

  // Different mean
  RooGaussian signalG1("signalG1","Signal PDF 1",JpsiMass,meanSig1,sigmaSig1);
  RooGaussian signalG2("signalG2","Signal PDF 2",JpsiMass,meanSig2,sigmaSig2);
  // Same mean
  RooGaussian signalG2OneMean("signalG2OneMean","Signal PDF 2",JpsiMass,meanSig1,sigmaSig2);

  RooRealVar coeffGauss("coeffGauss","Relative norm of the two signal gaussians",0.42,0.,1.);

  // Different mean
  RooAddPdf sigPDF("sigPDF","Total signal pdf",signalG1,signalG2,coeffGauss);
  // Same mean
  RooAddPdf sigPDFOneMean("sigPDFOneMean","Total signal pdf",signalG1,signalG2OneMean,coeffGauss);

  // SIGNAL: Crystal Ball shape
  RooRealVar alpha("alpha","#alpha of CB",0.96,0.,5.);
  RooRealVar enne("enne","n of CB",3.,0.,8.);

  RooRealVar NSig("NSig","Number of signal events",5000.,10.,10000000.);
  RooRealVar NBkg("NBkg","Number of background events",2000.,10.,10000000.);

  RooCBShape sigCB("sigCB","Signal CB PDF",JpsiMass,meanSig1,sigmaSig1,alpha,enne);

  RooAddPdf sigCBGauss("sigCBGauss","Signal CB+Gauss PDF",sigCB,signalG2,coeffGauss);

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
  // Total PDF (signal CB+Gauss)
  RooAddPdf totPDFCBGauss("totPDFCBGauss","Total pdf",RooArgList(sigCBGauss,CPolFunct2),RooArgList(NSig,NBkg));
  // Total PDF (signal CB+Gauss and exponential background)
  RooAddPdf totPDFExpCBGauss("totPDFExpCBGauss","Total pdf",RooArgList(sigCBGauss,expFunct),RooArgList(NSig,NBkg));

  //GG
  // fix some parameters 
  // alpha.setConstant(kTRUE); 
  // enne.setConstant(kTRUE); 

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  GGdata->setWeightVar(MCweight);
  RooDataSet *GGdataTr = (RooDataSet*)GGdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GGdataTr->setWeightVar(MCweight);

  if (sidebandPrefit) {
    PolFunct2.fitTo(*GGdata,Range("left,right"));
    CcoefPol1.setVal(coefPol1.getVal());
    CcoefPol1.setConstant(kTRUE);
    // float theCurrentVal = coefPol1.getVal();
    // coefPol1.setRange(theCurrentVal-0.01,theCurrentVal+0.01);
  }

  totPDFCBGauss.fitTo(*GGdata,Extended(1),Save(1),Minos(0));
  // totPDFExp.fitTo(*GGdata,Extended(1),Save(1));

  RooPlot *GGmframe = JpsiMass.frame();
  sprintf(reducestr,"Mass fit for glb-glb muons p_{T} = %s GeV,   #eta = %s",prange,etarange);
  GGmframe->SetTitle(reducestr);
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  // PolFunct2.plotOn(GGmframe,Range("left,right"));
  totPDFCBGauss.plotOn(GGmframe,Range("all",kTRUE));
  totPDFCBGauss.plotOn(GGmframe,Components(RooArgSet(CPolFunct2,sigCBGauss)),DrawOption("F"),FillColor(kGreen));
  totPDFCBGauss.plotOn(GGmframe,Components(CPolFunct2),DrawOption("F"),FillColor(kRed));
  // totPDFCBGaussExp.plotOn(GGmframe);
  // totPDFCBGaussExp.plotOn(GGmframe,Components(RooArgSet(expFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  // totPDFCBGaussExp.plotOn(GGmframe,Components(expFunct),DrawOption("F"),FillColor(kRed));
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  double NSigGG = NSig.getVal();   double errSigGG = NSig.getError();
  double resolGG = (coeffGauss.getVal()*coeffGauss.getVal()*sigmaSig1.getVal() + (1-coeffGauss.getVal())*(1-coeffGauss.getVal())*sigmaSig2.getVal()) / (coeffGauss.getVal()*coeffGauss.getVal() + (1-coeffGauss.getVal())*(1-coeffGauss.getVal())); 

  TCanvas c1;
  c1.cd();GGmframe->Draw();
  sprintf(reducestr,"GGmassfit_pT%s_eta%s.gif",prange,etarange);
  c1.SaveAs(reducestr);

  // fix some parameters 
  alpha.setConstant(kTRUE); 
  enne.setConstant(kTRUE); 
  // CcoefPol1.setVal(-1.0);
  // CcoefPol2.setConstant(kFALSE);

  //GT
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  GTdata->setWeightVar(MCweight);
  RooDataSet *GTdataTr = (RooDataSet*)GTdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GTdataTr->setWeightVar(MCweight);

  if (sidebandPrefit) {
    // coefPol1.setRange(-15000.,15000.);
    CcoefPol1.setConstant(kFALSE);
    PolFunct2.fitTo(*GTdata,Range("left,right"));
    CcoefPol1.setVal(coefPol1.getVal());
    CcoefPol1.setConstant(kTRUE);
    // float theCurrentVal = coefPol1.getVal();
    // coefPol1.setRange(theCurrentVal-0.01,theCurrentVal+0.01);
  }

  totPDFCBGauss.fitTo(*GTdata,Extended(1),Save(1),Minos(0));
  // totPDFExpCBGauss.fitTo(*GTdata,Extended(1),Save(1),Minos(0));

  RooPlot *GTmframe = JpsiMass.frame();
  sprintf(reducestr,"Mass fit for glb-trk muons p_{T} = %s GeV,  #eta = %s",prange,etarange);
  GTmframe->SetTitle(reducestr);
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  totPDFCBGauss.plotOn(GTmframe,Range("all",kTRUE));
  totPDFCBGauss.plotOn(GTmframe,Components(RooArgSet(CPolFunct2,sigCBGauss)),DrawOption("F"),FillColor(kGreen));
  totPDFCBGauss.plotOn(GTmframe,Components(CPolFunct2),DrawOption("F"),FillColor(kRed)); 
  /* totPDFExpCBGauss.plotOn(GTmframe,Range("all",kTRUE));
  totPDFExpCBGauss.plotOn(GTmframe,Components(RooArgSet(expFunct,sigCBGauss)),DrawOption("F"),FillColor(kGreen));
  totPDFExpCBGauss.plotOn(GTmframe,Components(expFunct),DrawOption("F"),FillColor(kRed)); */
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  double NSigGT = NSig.getVal();   double errSigGT = NSig.getError();
  double resolGT = (coeffGauss.getVal()*coeffGauss.getVal()*sigmaSig1.getVal() + (1-coeffGauss.getVal())*(1-coeffGauss.getVal())*sigmaSig2.getVal()) / (coeffGauss.getVal()*coeffGauss.getVal() + (1-coeffGauss.getVal())*(1-coeffGauss.getVal()));

  TCanvas c2;
  c2.cd();GTmframe->Draw();
  sprintf(reducestr,"GTmassfit_pT%s_eta%s.gif",prange,etarange);
  c2.SaveAs(reducestr); 

  //GC
  /* RooDataSet *GCdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GC");
  GCdata->setWeightVar(MCweight);
  RooDataSet *GCdataTr = (RooDataSet*)GCdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GCdataTr->setWeightVar(MCweight);
  
  meanSig1.setVal(3.1);
  meanSig1.setConstant(true);
  totPDFCBGauss.fitTo(*GCdata,Extended(1),Save(1),Minos(0));

  RooPlot *GCmframe = JpsiMass.frame();
  GCmframe->SetTitle("Mass fit for glb-trk muons");
  GCdata->plotOn(GCmframe,DataError(RooAbsData::SumW2));
  totPDFCBGauss.plotOn(GCmframe);
  totPDFCBGauss.plotOn(GCmframe,Components(RooArgSet(PolFunct2,sigPDFCBGauss)),DrawOption("F"),FillColor(kGreen));
  totPDFCBGauss.plotOn(GCmframe,Components(PolFunct2),DrawOption("F"),FillColor(kRed));
  GCdata->plotOn(GCmframe,DataError(RooAbsData::SumW2));
  double NSigGC = NSig.getVal();   double errSigGC = NSig.getError();

  TCanvas c3;
  c3.cd();GCmframe->Draw();
  c3.SaveAs("GCmassfit.gif"); */

  cout << endl << "pT = " << prange << " GeV; |eta| = " << etarange << endl;
  cout << endl << "GG J/psi yields:                 GT J/psi yields:" << endl;
  cout << "True MC : " << GGdataTr->numEntries(true) << "                   True MC : " << GTdataTr->numEntries(true) << endl;
  cout << "Fit : " << NSigGG << " +/- " << errSigGG << "        Fit : " << NSigGT << " +/- " << errSigGT << endl;
  cout << "Resolution : " << resolGG*1000. << " MeV         Resolution : " << resolGT*1000. << " MeV" << endl; 
  // cout << "GC J/psi yields:" << endl;
  // cout << "True MC : " << GCdataTr->numEntries(true) << " Fit : " << NSigGC << " +/- " << errSigGC << endl;

  return 1;
}
