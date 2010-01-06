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

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

void defineBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(JpsiMass,{CoefPol1[-0.05,-1500.,1500.],CcoefPol2[0.]})");

  //Exponential
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-5.,-9.,0.1])");

  return;
}

void defineSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  ws->factory("Gaussian::signalG1(JpsiMass,meanSig1[3.1,3.05,3.15],sigmaSig1[0.02,0.,0.2])");
  ws->factory("Gaussian::signalG2(JpsiMass,meanSig2[3.1,3.05,3.15],sigmaSig2[0.03,0.,0.2])");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(JpsiMass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(JpsiMass,meanSig1,sigmaSig1,alpha[0.96,0.,8.],enne[3.,0.,10.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.5,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  return;
}

void getrange(string &varRange, float *varmin, float *varmax)
{
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
   cout << varRange.c_str() << ": range not valid!" << endl;
    assert(0);
  }

 return;
}

void prefitSideband(RooWorkspace *ws, const int DataCat)
{
  ws->var("coefExp")->setConstant(kFALSE);

  RooDataSet *tmpdata;
  if(DataCat == 0) tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GG");
  else if(DataCat == 1) tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GT");
  else if(DataCat == 2) tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::TT");

  ws->pdf("expFunct")->fitTo(*tmpdata,Range("left,right"),SumW2Error(kTRUE));

  ws->var("coefExp")->setConstant(kTRUE);

  return;
}

void setRanges(RooWorkspace *ws)
{
  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;

  ws->var("JpsiMass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("JpsiMass")->setRange("left",JpsiMassMin,2.9);
  ws->var("JpsiMass")->setRange("right",3.3,JpsiMassMax);

  ws->cat("JpsiType")->setRange("glbglb","GG");
  ws->cat("JpsiType")->setRange("glbtrk","GT");
  ws->cat("JpsiType")->setRange("trktrk","TT");

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

void drawResults(RooWorkspace *ws, const int DataCat, const string prange, const string etarange)
{
  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");
  RooAbsPdf *background = ws->pdf("expFunct");

  char reducestr[200];

  RooPlot *mframe = ws->var("JpsiMass")->frame();

  if(DataCat == 0) sprintf(reducestr,"Mass fit for glb-glb muons p_{T} = %s GeV,   #eta = %s",prange.c_str(),etarange.c_str());
  else if(DataCat == 1) sprintf(reducestr,"Mass fit for glb-trk muons p_{T} = %s GeV,   #eta = %s",prange.c_str(),etarange.c_str());
  else if(DataCat == 2) sprintf(reducestr,"Mass fit for trk-trk muons p_{T} = %s GeV,   #eta = %s",prange.c_str(),etarange.c_str());

  mframe->SetTitle(reducestr);

  if(DataCat == 0) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else if(DataCat == 1) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  else if(DataCat == 2) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));

  totPDF->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components(*background),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));

  if(DataCat == 0) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else if(DataCat == 1) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  else if(DataCat == 2) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));

  TCanvas c1;
  c1.cd();mframe->Draw();
  if(DataCat == 0) sprintf(reducestr,"GGmassfit_pT%s_eta%s.gif",prange.c_str(),etarange.c_str());
  else if(DataCat == 1) sprintf(reducestr,"GTmassfit_pT%s_eta%s.gif",prange.c_str(),etarange.c_str());
  else if(DataCat == 2) sprintf(reducestr,"TTmassfit_pT%s_eta%s.gif",prange.c_str(),etarange.c_str());
  c1.SaveAs(reducestr);

  return;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol)
{
  Nsig   = ws->var("NSig")->getVal();
  errSig = ws->var("NSig")->getError();
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  const double sigmaSig2 = ws->var("sigmaSig2")->getVal();

  resol = (coeffGauss*coeffGauss*sigmaSig1 + (1-coeffGauss)*(1-coeffGauss)*sigmaSig2) / (coeffGauss*coeffGauss + (1-coeffGauss)*(1-coeffGauss)); 

  return;
}

int main(int argc, char* argv[])
{
  gROOT->SetStyle("Plain");

  char *filename;
  string prange;
  string etarange;
  bool sidebandPrefit = false;
  bool prefitSignalMass = false;

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

      case 'c':{
        prefitSignalMass = true;
        cout << "Signal MC pre-fitting activated" << endl;
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
  reddata->setWeightVar("MCweight");

  ws->import(*reddata);

  setRanges(ws);

  //DEFINE SIGNAL AND BACKGROUND
  defineSignal(ws);
  defineBackground(ws);

  // Total PDF (signal CB+Gauss)
  ws->factory("SUM::totPDF(NSig[5000.,10.,10000000.]*sigCBGauss,NBkg[2000.,10.,10000000.]*expFunct)");

  //Make subsamples to be used later

  RooDataSet *GGdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *GTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *TTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT && (MCType == MCType::PR || MCType == MCType::NP)");

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  RooDataSet *TTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT");

  //GG CASE

  // fix some parameters 
  // ws->var("alpha")->setConstant(kTRUE); 
  // ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,0);

  if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*GGdataTr,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);
  }

  ws->pdf("totPDF")->fitTo(*GGdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,0,prange,etarange);

  double NSigGG, errSigGG,resolGG;
  printResults(ws,NSigGG,errSigGG,resolGG);

  //GT case

  // fix some parameters 
  //ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,1);

  ws->pdf("totPDF")->fitTo(*GTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,1,prange,etarange);

  double NSigGT, errSigGT,resolGT;
  printResults(ws,NSigGT,errSigGT,resolGT);

  //TT case

  // fix some parameters 
  //ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,2);

  ws->pdf("totPDF")->fitTo(*TTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,2,prange,etarange);

  double NSigTT, errSigTT,resolTT;
  printResults(ws,NSigTT,errSigTT,resolTT);

  cout << endl << "pT = " << prange << " GeV; |eta| = " << etarange << endl;
  cout << endl << "GG J/psi yields:                 GT J/psi yields:                   TT J/psi yields:" << endl;
  cout << "True MC : " << GGdataTr->sumEntries() << "                   True MC : " << GTdataTr->sumEntries() << "                   True MC : " << TTdataTr->sumEntries() << endl;
  cout << "Fit : " << NSigGG << " +/- " << errSigGG << "        Fit : " << NSigGT << " +/- " << errSigGT << "        Fit : " << NSigTT << " +/- " << errSigTT << endl;
  cout << "Resolution : " << resolGG*1000. << " MeV         Resolution : " << resolGT*1000. << " MeV        Resolution : " << resolTT*1000. << " MeV" << endl; 

  return 1;
}
