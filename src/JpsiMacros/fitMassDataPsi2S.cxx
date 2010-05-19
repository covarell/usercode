// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
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
  ws->factory("Gaussian::signalG1(JpsiMass,meanSig1[3.1,3.05,3.15],sigmaSig1[0.07,0.027,0.2])");
  ws->factory("Gaussian::signalG2(JpsiMass,meanSig2[3.7,3.6,3.8],sigmaSig1)");
  // ws->factory("Gaussian::signalG2(JpsiMass,meanSig1,sigmaSig1");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(JpsiMass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(JpsiMass,meanSig1,sigmaSig1,alpha[0.5,0.,2.],enne[10.,2.,50.])");

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
  else if(DataCat == 3) tmpdata = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GT || JpsiType == JpsiType::GG");

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

  /* ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");*/

  return;
}

void drawResults(RooWorkspace *ws, const int DataCat)
{

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");

  char reducestr[200];

  ws->var("JpsiMass")->SetTitle("#mu^{+} #mu^{-} mass");
  RooPlot *mframe = ws->var("JpsiMass")->frame();

  if(DataCat == 0) sprintf(reducestr,"Mass fit for glb-glb muons");
  else if(DataCat == 1) sprintf(reducestr,"Mass fit for glb-trk muons");
  else if(DataCat == 2) sprintf(reducestr,"Mass fit for trk-trk muons");
  else if(DataCat == 3) sprintf(reducestr,"Mass fit for all muon types");

  mframe->SetTitle(reducestr);

  if(DataCat == 0) {
    /* data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"),MarkerColor(4));
       data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));*/
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"),Binning(32));
  }
  else if(DataCat == 1) {
    /* data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"),MarkerColor(4));
       data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));*/
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  }
  else if(DataCat == 2) { 
    /* data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && MCType == MCType::BK"),MarkerColor(4));
       data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));*/
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));
  }
  else if(DataCat == 3) { 
    // data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG || JpsiType == JpsiType::GT"),Binning(20));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(32));
  }

  totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

  // mframe->addObject(pl1);

  // if(DataCat == 0) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  // else if(DataCat == 1) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  // else if(DataCat == 2) data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));

  TCanvas c1;
  c1.cd(); /* c1.SetLogy(1) */; mframe->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  // t->SetTextFont(63);
  t->SetTextSizePixels(22);
  t->DrawLatex(0.6,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");

  t->SetTextAlign(13); //align at top left
  t->SetTextAlign(12); // left, vertically centered
  t->SetTextAlign(22); // centered horizontally and vertically
  t->SetTextAlign(11); //default bottom alignment

  if(DataCat == 0) sprintf(reducestr,"GGmassfit.gif");
  else if(DataCat == 1) sprintf(reducestr,"GTmassfit.gif");
  else if(DataCat == 2) sprintf(reducestr,"TTmassfit.gif");
  // else if(DataCat == 3) sprintf(reducestr,"GGTmassfit.gif");
  else if(DataCat == 3) sprintf(reducestr,"allmassfit.gif");
  c1.SaveAs(reducestr);

  TCanvas c2;
  ws->var("JpsiPt")->setRange(0.,20.);
  ws->var("JpsiPt")->SetTitle("J/#psi p_{T}");
  RooPlot *ptframe = ws->var("JpsiPt")->frame();
  data->plotOn(ptframe,DataError(RooAbsData::SumW2),Cut("JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20),LineColor(kBlue),MarkerColor(kBlue));
  data->plotOn(ptframe,DataError(RooAbsData::SumW2),Cut("(JpsiType == JpsiType::GG || JpsiType == JpsiType::GT) && JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20));
  data->plotOn(ptframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20),LineColor(kRed),MarkerColor(kRed));
  c2.cd(); /* c1.SetLogy(1) */; ptframe->Draw();
  // t->DrawLatex(0.3,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");
  c2.SaveAs("allptplot.gif");

  TCanvas c3;
  ws->var("JpsiEta")->SetTitle("J/#psi y");
  RooPlot *etaframe = ws->var("JpsiEta")->frame();
  data->plotOn(etaframe,DataError(RooAbsData::SumW2),Cut("JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20),LineColor(kBlue),MarkerColor(kBlue));
  data->plotOn(etaframe,DataError(RooAbsData::SumW2),Cut("(JpsiType == JpsiType::GG || JpsiType == JpsiType::GT) && JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20));
  data->plotOn(etaframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(20),LineColor(kRed),MarkerColor(kRed));
  c3.cd(); /* c1.SetLogy(1) */; etaframe->Draw();
  // t->DrawLatex(0.2,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");
  c3.SaveAs("alletaplot.gif");
  
  TCanvas c4;
  RooPlot *ctframe = ws->var("Jpsict")->frame();
  // data->plotOn(ctframe,DataError(RooAbsData::SumW2),Cut("(JpsiType == JpsiType::GG || JpsiType == JpsiType::GT) && JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(18));
  data->plotOn(ctframe,DataError(RooAbsData::SumW2),Cut("JpsiMass < 3.8 && JpsiMass > 3.6"),Binning(18));
  c4.cd();  c4.SetLogy(1) ; ctframe->Draw();
  c4.SaveAs("allctauplot.gif");

  return;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol)
{
  Nsig   = ws->var("NSig")->getVal();
  errSig = ws->var("NSig")->getError();
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  // const double sigmaSig2 = ws->var("sigmaSig2")->getVal();

  // resol = (coeffGauss*coeffGauss*sigmaSig1 + (1-coeffGauss)*(1-coeffGauss)*sigmaSig2) / (coeffGauss*coeffGauss + (1-coeffGauss)*(1-coeffGauss)); 
  
  resol = sigmaSig1;

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
  // reddata->setWeightVar("MCweight");

  ws->import(*reddata);

  setRanges(ws);

  //DEFINE SIGNAL AND BACKGROUND
  defineSignal(ws);
  defineBackground(ws);

  // Total PDF (signal CB+Gauss)
  ws->factory("SUM::totPDF(NSig[5000.,4.,10000000.]*signalG1,NSig2[5000.,-10.,10000000.]*signalG2,NBkg[2000.,4.,10000000.]*expFunct)");

  //Make subsamples to be used later
  /*
  RooDataSet *GGdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *GTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *TTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT && (MCType == MCType::PR || MCType == MCType::NP)");*/

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  RooDataSet *GGTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG || JpsiType == JpsiType::GT");
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  RooDataSet *TTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT");
  /* ws->var("JpsiMass")->setBins(120);

  RooDataHist *GGdataBin = new RooDataHist("GGdataBin","GGdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdata);
  RooDataHist *GTdataBin = new RooDataHist("GTdataBin","GTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GTdata);
  RooDataHist *TTdataBin = new RooDataHist("TTdataBin","TTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*TTdata);*/

  //GG CASE

  // fix some parameters 
  ws->var("alpha")->setConstant(kTRUE); 
  ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,3);

  /* if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*GGdataTr,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);
    }*/

  // ws->var("coeffGauss")->setVal(1.0);
  // ws->var("coeffGauss")->setConstant(kTRUE);

  ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  // ws->pdf("totPDF")->fitTo(*GGdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  // ws->pdf("totPDF")->fitTo(*GGdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,3);
  // drawResults(ws,0);

  double NSigGG, errSigGG,resolGG;
  printResults(ws,NSigGG,errSigGG,resolGG);

  //GT case

  // fix some parameters 
  /* ws->var("alpha")->setConstant(kTRUE); 
  ws->var("enne")->setConstant(kTRUE); 
  // ws->var("meanSig1")->setConstant(kTRUE);

  if (sidebandPrefit) prefitSideband(ws,1);

  // ws->pdf("totPDF")->fitTo(*GGTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF")->fitTo(*GTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,1);

  double NSigGT, errSigGT,resolGT;
  printResults(ws,NSigGT,errSigGT,resolGT);

  //TT case

  // fix some parameters 
  //ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,2);

  // ws->pdf("totPDF")->fitTo(*TTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF")->fitTo(*TTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,2);

  double NSigTT, errSigTT,resolTT;
  printResults(ws,NSigTT,errSigTT,resolTT);*/

  cout << endl << "GG J/psi yields:                 GT J/psi yields:                   TT J/psi yields:" << endl;
  // cout << "True MC : " << GGdataTr->sumEntries() << "                   True MC : " << GTdataTr->sumEntries() << "                   True MC : " << TTdataTr->sumEntries() << endl;
  // cout << "Fit : " << NSigGG << " +/- " << errSigGG << "        Fit : " << NSigGT << " +/- " << errSigGT << "        Fit : " << NSigTT << " +/- " << errSigTT << endl;
  cout << "Fit : " << NSigGG << " +/- " << errSigGG << endl;
  cout << "Resolution : " << resolGG*1000. << " MeV" << endl; 

  return 1;
}
