// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TLatex.h>
#include <TMath.h>
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
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooWorkspace.h"

using namespace RooFit;
TCanvas* c1;

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

void defineBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(JpsiMass,{CoefPol1[-0.05,-1500.,1500.],CcoefPol2[0.]})");

  //Exponential
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-1.,-3.,0.1])");

  return;
}

void defineSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  ws->factory("Gaussian::signalG1(JpsiMass,meanSig1[3.1,3.05,3.15],sigmaSig1[0.02,0.008,0.2])");
  ws->factory("Gaussian::signalG2(JpsiMass,meanSig2[3.1,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(JpsiMass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(JpsiMass,meanSig1,sigmaSig1,alpha[0.5,0.,3.],enne[10.,1.,30.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.5,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGaussOneMean(coeffGauss*sigCB,signalG1)");

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
  else if(DataCat == 3) tmpdata = (RooDataSet*)ws->data("data");

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

void drawResults(RooWorkspace *ws, const int DataCat, const string prange, const string etarange, const int nFitPar)
{

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");

  char reducestr[200];

  ws->var("JpsiMass")->SetTitle("#mu^{+} #mu^{-} mass");
  RooPlot *mframe = ws->var("JpsiMass")->frame();

  if(DataCat == 0) sprintf(reducestr,"Mass fit for glb-glb muons p_{T} = %s GeV,   |y| = %s",prange.c_str(),etarange.c_str()); 
  else if(DataCat == 1) sprintf(reducestr,"Mass fit for glb-trk muons p_{T} = %s GeV,   |y| = %s",prange.c_str(),etarange.c_str());
  else if(DataCat == 2) sprintf(reducestr,"Mass fit for trk-trk muons p_{T} = %s GeV,   |y| = %s",prange.c_str(),etarange.c_str());
  else if(DataCat == 3) sprintf(reducestr,"Mass fit for all muons p_{T} = %s GeV,   |y| = %s",prange.c_str(),etarange.c_str());

  mframe->SetTitle(reducestr);
  RooHist* hresid;
  double chi2;

  if(DataCat == 0) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = mframe->residHist();
    hresid->SetName("hresid");
    chi2 = mframe->chiSquare(nFitPar);
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  }
  else if(DataCat == 1) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = mframe->residHist();
    hresid->SetName("hresid");
    chi2 = mframe->chiSquare(nFitPar);
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));
  }
  else if(DataCat == 2) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = mframe->residHist();
    hresid->SetName("hresid");
    chi2 = mframe->chiSquare(nFitPar);
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT && (MCType == MCType::BK || MCType == MCType::NP)"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::TT"));
  }
  else if(DataCat == 3) {
    data->plotOn(mframe,DataError(RooAbsData::SumW2));
    totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = mframe->residHist();
    hresid->SetName("hresid");
    chi2 = mframe->chiSquare(nFitPar);
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::BK"),MarkerColor(4));
    data->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::BK || MCType == MCType::NP"),MarkerColor(2));
    data->plotOn(mframe,DataError(RooAbsData::SumW2));
  }

  totPDF->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

  // NO RESIDUALS 
  c1 = new TCanvas();
  c1->cd();  /* c1.SetLogy(1) */; 
  mframe->Draw();
  // Do not use ugly paramOn layout, but nice TLatex
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextAlign(12); // left, vertically centered
  // t->SetTextFont(63);
  t->SetTextSize(0.036);
  sprintf(reducestr,"N_{sig} = %.0f #pm %.0f",int(ws->var("NSig")->getVal()/10)*10.,int(ws->var("NSig")->getError()/10)*10.);
  t->DrawLatex(0.15,0.9,reducestr);
  // sprintf(reducestr,"#LT m #GT = %.1f #pm %.1f MeV/c^{2}",ws->var("meanSig1")->getVal()*1000.,ws->var("meanSig1")->getError()*1000.);
  sprintf(reducestr,"#LT m #GT = 3.097.1 #pm 0.3 MeV/c^{2}");
  t->DrawLatex(0.15,0.86,reducestr);
  sprintf(reducestr,"#sigma_{1} = %.1f #pm %.1f MeV/c^{2}",ws->var("sigmaSig1")->getVal()*1000.,ws->var("sigmaSig1")->getError()*1000.);
  t->DrawLatex(0.15,0.82,reducestr);
  sprintf(reducestr,"#sigma_{2} = %.1f #pm %.1f MeV/c^{2}",ws->var("sigmaSig2")->getVal()*1000.,ws->var("sigmaSig2")->getError()*1000.);
  t->DrawLatex(0.15,0.78,reducestr);

  // t->SetTextAlign(13); //align at top left
  t->SetTextAlign(12); // left, vertically centered
  //  t->SetTextAlign(22); // centered horizontally and vertically
  //  t->SetTextAlign(11); //default bottom alignment
  
  // WITH RESIDUALS 
  /* c1 = new TCanvas("c1","The Canvas",200,10,600,880);
  c1->cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.30);
  pad2->Draw();

  pad1->cd(); mframe->Draw();

  RooPlot* mframe2 =  ws->var("JpsiMass")->frame(Title("Residuals Distribution")) ;
  mframe2->addPlotable(hresid,"P") ;  

  pad2->cd(); mframe2->Draw();

  int nDOF = ws->var("JpsiMass")->getBinning().numBins() - nFitPar;

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  // t->SetTextFont(63);
  t->SetTextSizePixels(22);
  sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  if (chi2 < 10.) t->DrawLatex(0.6,0.9,reducestr);

  t->SetTextAlign(13); //align at top left
  t->SetTextAlign(12); // left, vertically centered
  t->SetTextAlign(22); // centered horizontally and vertically
  t->SetTextAlign(11); //default bottom alignment
  */
  
  c1->Update();

  if(DataCat == 0) sprintf(reducestr,"pictures/GGmassfit_pT%s_eta%s.pdf",prange.c_str(),etarange.c_str());
  else if(DataCat == 1) sprintf(reducestr,"pictures/GTmassfit_pT%s_eta%s.pdf",prange.c_str(),etarange.c_str());
  else if(DataCat == 2) sprintf(reducestr,"pictures/TTmassfit_pT%s_eta%s.pdf",prange.c_str(),etarange.c_str());
  else if(DataCat == 3) sprintf(reducestr,"pictures/ALLmassfit_pT%s_eta%s.pdf",prange.c_str(),etarange.c_str());

  c1->SaveAs(reducestr);

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
  errresol = (0.5/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));

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
  ws->var("JpsiMass")->setBins(60);

  RooDataSet *reddataTr = (RooDataSet*)reddata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  RooDataSet *GGdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *GTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT && (MCType == MCType::PR || MCType == MCType::NP)");
  RooDataSet *TTdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT && (MCType == MCType::PR || MCType == MCType::NP)");

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  RooDataSet *TTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT");

  RooDataHist *reddataBin = new RooDataHist("reddataBin","reddataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*reddata);

  RooDataHist *GGdataBin = new RooDataHist("GGdataBin","GGdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdata);

  cout << "Number of events to fit  = " << GGdata->sumEntries() << endl; 

  RooDataHist *GTdataBin = new RooDataHist("GTdataBin","GTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GTdata);
  RooDataHist *TTdataBin = new RooDataHist("TTdataBin","TTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*TTdata);

  RooDataHist *GGdataTrBin = new RooDataHist("GGdataTrBin","GGdataTrBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdataTr);

  cout << "Number of true events to fit  = " << GGdataTr->sumEntries() << endl;

  // OPTION A: All together

  /* ws->var("JpsiMass")->setBins(90);
  // if (etamin < 1.0)  ws->var("JpsiMass")->setBins(18);

  // fix some parameters 
  // ws->var("alpha")->setConstant(kTRUE); 
  // ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,0);

  if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*GGdataTrBin,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    // ws->var("alpha")->setConstant(kTRUE);
  }

  // ws->var("coeffGauss")->setVal(1.0);
  // ws->var("coeffGauss")->setConstant(kTRUE);

  // ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  RooFitResult *rfr = ws->pdf("totPDF")->fitTo(*reddataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,3,prange,etarange,rfr->floatParsFinal().getSize());

  double NSig, errSig,resol,errresol;
  printResults(ws,NSig,errSig,resol,errresol);

  char oFile[200];
  sprintf(oFile,"results/results_pT%s_eta%s.txt",prange.c_str(),etarange.c_str());
  ofstream outputFile(oFile);
  outputFile << "AL " << reddataTr->sumEntries() << " " << NSig << " " << errSig << endl;
  outputFile << "RE " << 0. << " " << resol*1000. << " " << errresol*1000. << endl;
  outputFile << endl;*/

  // OPTION B: Separate

  //GG CASE
  // if (etamin < 1.0)  ws->var("JpsiMass")->setBins(18);

  // fix some parameters 
  // ws->var("alpha")->setConstant(kTRUE); 
  // ws->var("enne")->setConstant(kTRUE); 

  if (sidebandPrefit) prefitSideband(ws,0);

  if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*GGdataTrBin,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    // ws->var("alpha")->setConstant(kTRUE);
  }

  // ws->var("coeffGauss")->setVal(1.0);
  // ws->var("coeffGauss")->setConstant(kTRUE);

  // ws->pdf("totPDF")->fitTo(*GGdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  RooFitResult *rfr = ws->pdf("totPDF")->fitTo(*GGdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,0,prange,etarange,rfr->floatParsFinal().getSize());

  double NSigGG, errSigGG,resolGG,errresolGG;
  printResults(ws,NSigGG,errSigGG,resolGG,errresolGG);

  cout << "GG " << GGdataTr->sumEntries() << " " << NSigGG << " " << errSigGG << endl;

  /* 
  //GT case
  // ws->var("JpsiMass")->setBins(45);

  // fix some parameters for barrel
  if (etamin < 1.0) {
    ws->var("alpha")->setConstant(kTRUE); 
    ws->var("enne")->setConstant(kTRUE);
  } 
  // ws->var("meanSig1")->setConstant(kTRUE);

  if (sidebandPrefit) prefitSideband(ws,1);

  // ws->pdf("totPDF")->fitTo(*GTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  rfr = ws->pdf("totPDF")->fitTo(*GTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,1,prange,etarange,rfr->floatParsFinal().getSize());

  double NSigGT, errSigGT,resolGT,errresolGT;
  printResults(ws,NSigGT,errSigGT,resolGT,errresolGT);

  //TT case

  // fix some parameters 
  ws->var("alpha")->setConstant(kTRUE); 
  ws->var("enne")->setConstant(kTRUE);
  ws->var("meanSig1")->setConstant(kTRUE);
  ws->var("sigmaSig1")->setConstant(kTRUE);

  if (sidebandPrefit) prefitSideband(ws,2);

  // ws->pdf("totPDF")->fitTo(*TTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  // rfr = ws->pdf("totPDF")->fitTo(*TTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  // drawResults(ws,2,prange,etarange,rfr->floatParsFinal().getSize());

  // double NSigTT, errSigTT,resolTT,errresolTT;
  // printResults(ws,NSigTT,errSigTT,resolTT,errresolTT);

  char oFile[200];
  sprintf(oFile,"results/CBGaussCommonMean/results_pT%s_eta%s.txt",prange.c_str(),etarange.c_str());
  ofstream outputFile(oFile);
  outputFile << "GG " << GGdataTr->sumEntries() << " " << NSigGG << " " << errSigGG << endl;
  outputFile << "RE " << 0. << " " << resolGG*1000. << " " << errresolGG*1000. << endl;
  outputFile << "GT " << GTdataTr->sumEntries() << " " << NSigGT << " " << errSigGT << endl;
  // outputFile << "TT " << TTdataTr->sumEntries() << " " << NSigTT << " " << errSigTT << endl;
  outputFile << endl; */

  return 1;
}
