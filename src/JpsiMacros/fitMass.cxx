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
  JpsiMass.setRange("left",2.6,2.9);
  JpsiMass.setRange("right",3.3,3.5);

  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi pt",0.,200.,"GeV/c");
  RooRealVar* JpsiEta = new RooRealVar("JpsiEta","J/psi eta",-2.7,2.7);
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,5.);

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
  RooRealVar coefPol1("coefPol1","linear coefficient of bkg PDF",0.,-9.,9.);
  RooRealVar coefPol2("coefPol2","quadratic coefficient of bkg PDF",0.,-1.,1.);
  RooPolynomial GGPolFunct("GGPolFunc","GGPolFunc",JpsiMass,coefPol1);
  RooPolynomial GCPolFunct("GCPolFunc","GCPolFunc",JpsiMass,RooArgList(coefPol1,coefPol2));
 
  // BKG : exponential
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.);
  RooExponential GTexpFunct("GTexpFunc","GTexpFunc",JpsiMass,coefExp); 

  // SIGNAL: 2-Gaussians
  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.1,3.05,3.15);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.02,0.,0.5);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.1,3.05,3.15);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.04,0.,0.5);

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

  RooRealVar NSig("NSig","Number of signal events",5000.,10.,10000000.);
  RooRealVar NBkg("NBkg","Number of background events",1300.,10.,10000000.);

  // Total PDF
  RooAddPdf totPDF("totPDF","Total pdf",RooArgList(sigPDF,GGPolFunct),RooArgList(NSig,NBkg));
  // Total PDF (exponential background)
  RooAddPdf totPDFExp("totPDFexp","Total pdf",RooArgList(sigPDF,GTexpFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean)
  RooAddPdf totPDFOneMean("totPDFOneMean","Total pdf",RooArgList(sigPDFOneMean,GGPolFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean and exponential background)
  RooAddPdf totPDFExpOneMean("totPDFexpOneMean","Total pdf",RooArgList(sigPDFOneMean,GTexpFunct),RooArgList(NSig,NBkg));
  // Total PDF (signal Gaussians with one mean and quadratic background)
  RooAddPdf totPDF2Pol("totPDF2Pol","Total pdf",RooArgList(sigPDFOneMean,GCPolFunct),RooArgList(NSig,NBkg));

  //GG
  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  GGdata->setWeightVar(MCweight);
  RooDataSet *GGdataTr = (RooDataSet*)GGdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GGdataTr->setWeightVar(MCweight);

  totPDF.fitTo(*GGdata,Extended(1),Save(1));
  // totPDFExp.fitTo(*GGdata,Extended(1),Save(1));

  RooPlot *GGmframe = JpsiMass.frame();
  GGmframe->SetTitle("Mass fit for glb-glb muons");
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GGmframe);
  totPDF.plotOn(GGmframe,Components(RooArgSet(GGPolFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  totPDF.plotOn(GGmframe,Components(GGPolFunct),DrawOption("F"),FillColor(kRed));
  // totPDFExp.plotOn(GGmframe);
  // totPDFExp.plotOn(GGmframe,Components(RooArgSet(GTexpFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  // totPDFExp.plotOn(GGmframe,Components(GTexpFunct),DrawOption("F"),FillColor(kRed));
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  double NSigGG = NSig.getVal();   double errSigGG = NSig.getError();

  TCanvas c1;
  c1.cd();GGmframe->Draw();
  sprintf(reducestr,"GGmassfit_pT%s_eta%s.gif",prange,etarange);
  c1.SaveAs(reducestr);

  //GT
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  GTdata->setWeightVar(MCweight);
  RooDataSet *GTdataTr = (RooDataSet*)GTdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GTdataTr->setWeightVar(MCweight);

  totPDFExpOneMean.fitTo(*GTdata,Extended(1),Save(1)/* ,Minos(1)*/);

  RooPlot *GTmframe = JpsiMass.frame();
  GTmframe->SetTitle("Mass fit for glb-trk muons");
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  totPDFExpOneMean.plotOn(GTmframe);
  totPDFExpOneMean.plotOn(GTmframe,Components(RooArgSet(GTexpFunct,sigPDFOneMean)),DrawOption("F"),FillColor(kGreen));
  totPDFExpOneMean.plotOn(GTmframe,Components(GTexpFunct),DrawOption("F"),FillColor(kRed));
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  double NSigGT = NSig.getVal();   double errSigGT = NSig.getError();

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
  totPDFOneMean.fitTo(*GCdata,Extended(1),Save(1),Minos(0));

  RooPlot *GCmframe = JpsiMass.frame();
  GCmframe->SetTitle("Mass fit for glb-trk muons");
  GCdata->plotOn(GCmframe,DataError(RooAbsData::SumW2));
  totPDFOneMean.plotOn(GCmframe);
  totPDFOneMean.plotOn(GCmframe,Components(RooArgSet(GGPolFunct,sigPDFOneMean)),DrawOption("F"),FillColor(kGreen));
  totPDFOneMean.plotOn(GCmframe,Components(GGPolFunct),DrawOption("F"),FillColor(kRed));
  GCdata->plotOn(GCmframe,DataError(RooAbsData::SumW2));
  double NSigGC = NSig.getVal();   double errSigGC = NSig.getError();

  TCanvas c3;
  c3.cd();GCmframe->Draw();
  c3.SaveAs("GCmassfit.gif"); */

  cout << endl << "GG J/psi yields:" << endl;
  cout << "True MC : " << GGdataTr->numEntries(true) << " Fit : " << NSigGG << " +/- " << errSigGG << endl;
  cout << "GT J/psi yields:" << endl;
  cout << "True MC : " << GTdataTr->numEntries(true) << " Fit : " << NSigGT << " +/- " << errSigGT << endl;
  // cout << "GC J/psi yields:" << endl;
  // cout << "True MC : " << GCdataTr->numEntries(true) << " Fit : " << NSigGC << " +/- " << errSigGC << endl;

  return 1;
}
