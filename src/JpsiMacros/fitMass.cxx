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
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCurve.h"
#include "RooCategory.h"
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
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-1.,-2.,0.1])");

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
  ws->factory("CBShape::sigCB(JpsiMass,meanSig1,sigmaSig1,alpha[0.5,0.,2.],enne[20.,2.,50.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.5,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

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

void drawResults(RooWorkspace *ws, const int DataCat)
{
  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");

  char reducestr[200];

  ws->var("JpsiMass")->SetTitle("#mu^{+} #mu^{-} mass");
  RooPlot *mframe = ws->var("JpsiMass")->frame();

  if(DataCat == 0) sprintf(reducestr,"Mass fit for glb-glb muons");
  else if(DataCat == 1) sprintf(reducestr,"Mass fit for glb-trk muons");
  else if(DataCat == 2) sprintf(reducestr,"Mass fit for trk-trk muons");

  mframe->SetTitle(reducestr);
  RooHist* hresid;

  if(DataCat == 0) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = mframe->pullHist();
    hresid->SetName("hresid");
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  }
  else if(DataCat == 1) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  }
  else if(DataCat == 2) { 
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));
  }

  totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c1;
  TCanvas c2("c2","c2",10,10,700,300);
  if (DataCat == 0) {
    c2.cd();
    RooPlot* mframe2 =  ws->var("JpsiMass")->frame(Title("Residual Distribution")) ;
    mframe2->addPlotable(hresid,"P") ;
    mframe2->Draw();
    c2.SaveAs("testpull.gif");
  }
  c1.cd(); /* c1.SetLogy(1) */; mframe->Draw();
  if(DataCat == 0) sprintf(reducestr,"GGmassfit.gif");
  else if(DataCat == 1) sprintf(reducestr,"GTmassfit.gif");
  else if(DataCat == 2) sprintf(reducestr,"TTmassfit.gif");
  c1.SaveAs(reducestr);

  return;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol,double &errresol)
{
  Nsig   = ws->var("NSig")->getVal();
  errSig = ws->var("NSig")->getError();
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  const double sigmaSig2 = ws->var("sigmaSig2")->getVal();
  const double ecoeffGauss = ws->var("coeffGauss")->getError();
  const double esigmaSig1 = ws->var("sigmaSig1")->getError();
  const double esigmaSig2 = ws->var("sigmaSig2")->getError();

  resol = sqrt(coeffGauss*sigmaSig1*sigmaSig1 + (1-coeffGauss)*sigmaSig2*sigmaSig2);
  errresol = (1/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));

  return;
}

int main(int argc, char* argv[])
{
  gROOT->SetStyle("Plain");

  char *filename;
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

  RooDataSet *reddata = (RooDataSet*)fIn.Get("data");
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
  ws->var("JpsiMass")->setBins(120);

  RooDataHist *GGdataBin = new RooDataHist("GGdataBin","GGdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdata);
  RooDataHist *GTdataBin = new RooDataHist("GTdataBin","GTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GTdata);
  RooDataHist *TTdataBin = new RooDataHist("TTdataBin","TTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*TTdata);

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

  // ws->var("coeffGauss")->setVal(1.0);
  // ws->var("coeffGauss")->setConstant(kTRUE);

  // ws->pdf("totPDF")->fitTo(*GGdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF")->fitTo(*GGdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,0);

  double NSigGG, errSigGG,resolGG,errresolGG;
  printResults(ws,NSigGG,errSigGG,resolGG,errresolGG);

  //GT case

  // fix some parameters 
  ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 
  ws->var("meanSig1")->setConstant(kTRUE);

  if (sidebandPrefit) prefitSideband(ws,1);

  // ws->pdf("totPDF")->fitTo(*GTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF")->fitTo(*GTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,1);

  double NSigGT, errSigGT,resolGT,errresolGT;
  printResults(ws,NSigGT,errSigGT,resolGT,errresolGT);

  //TT case

  // fix some parameters 
  //ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,2);

  // ws->pdf("totPDF")->fitTo(*TTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF")->fitTo(*TTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,2);

  double NSigTT, errSigTT,resolTT,errresolTT;
  printResults(ws,NSigTT,errSigTT,resolTT,errresolTT);

  cout << endl << "GG J/psi yields:                 GT J/psi yields:                   TT J/psi yields:" << endl;
  cout << "True MC : " << GGdataTr->sumEntries() << "                   True MC : " << GTdataTr->sumEntries() << "                   True MC : " << TTdataTr->sumEntries() << endl;
  cout << "Fit : " << NSigGG << " +/- " << errSigGG << "        Fit : " << NSigGT << " +/- " << errSigGT << "        Fit : " << NSigTT << " +/- " << errSigTT << endl;
  cout << "Resolution : " << resolGG*1000. << " +/- " << errresolGG*1000. << " MeV  Resolution : " << resolGT*1000. << " +/- " << errresolGT*1000. << " MeV   Resolution : " << resolTT*1000. << " +/- " << errresolTT*1000. << " MeV" << endl; 

  return 1;
}
