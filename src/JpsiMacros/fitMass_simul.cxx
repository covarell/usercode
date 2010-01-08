// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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

static const Int_t nbinspt = 11;
static const Int_t nbinseta = 2;

double ptbinlimits[nbinspt+1];
double etabinlimits[nbinseta+1];

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

string nameCategory(string nameBase, int i)
{
  ostringstream ost;
  ost << nameBase << i;
  
  return ost.str();
}

string nameSplitVar(string nameBaseVar, int ip = -1, int ieta = -1)
{
  ostringstream ostP;
  ostP << "P" << ip;
  ostringstream ostEta;
  ostEta << "E" << ieta;

  string catstring_pt = ostP.str();
  string catstring_eta = ostEta.str(); 
  string result;

  if (ip == -1) result = nameBaseVar + "_" + catstring_eta;
  else if (ieta == -1) result = nameBaseVar + "_" + catstring_pt;
  else result = nameBaseVar + "_{" + catstring_pt + ";" + catstring_eta + "}";
  
  return result;
}

string nameCategoryCut(int ip = -1, int ieta = -1)
{
  ostringstream ostP2;
  ostP2 << "P" << ip;
  ostringstream ostEta2;
  ostEta2 << "E" << ieta;

  string catstring_pt = ostP2.str();
  string catstring_eta = ostEta2.str(); 
  string result;

  if (ip == -1) result = "JpsiEtaType == JpsiEtaType::" + catstring_eta;
  else if (ieta == -1) result = "JpsiPtType == JpsiPtType::" + catstring_pt;
  else result = "JpsiEtaType == JpsiEtaType::" + catstring_eta + " && JpsiPtType == JpsiPtType::" + catstring_pt;

  return result;
}

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

void prefitSideband(RooWorkspace *ws, const int DataCat)
{

  for(int i = 1;i<=nbinspt;i++){
    for(int j = 1;j<=nbinseta;j++){

      string myName = nameSplitVar("coefExp",i,j);
      if (ws->var(myName.c_str())) ws->var(myName.c_str())->setConstant(kFALSE);
     
      string halfcutstring = nameCategoryCut(i,j);
      string cutstring;

      if(DataCat == 0) cutstring = halfcutstring + " && JpsiType == JpsiType::GG";
      else if(DataCat == 1) cutstring = halfcutstring + " && JpsiType == JpsiType::GT"; 
      else if(DataCat == 2)  cutstring = halfcutstring + " && JpsiType == JpsiType::GG";

      RooDataSet *tmpdata = (RooDataSet*)ws->data("data")->reduce(cutstring.c_str());

      RooDataHist *tmpdataBin = new RooDataHist("tmpdataBin","tmpdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*tmpdata);
      
      cutstring = nameSplitVar("expFunct",i,j);
      ws->pdf(cutstring.c_str())->fitTo(*tmpdata,Range("left,right"),SumW2Error(kTRUE));
      // ws->pdf(cutstring.c_str())->fitTo(*tmpdataBin,Range("left,right"),SumW2Error(kTRUE));
  
      if (ws->var(myName.c_str())) ws->var(myName.c_str())->setConstant(kTRUE);
   
    }
  }

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
  RooAbsPdf *totalPDF = ws->pdf("totPDF_sim");
  RooCategory *JpsiPtType = ws->cat("JpsiPtType");
  RooCategory *JpsiEtaType = ws->cat("JpsiEtaType");
  RooRealVar *JpsiMass = ws->var("JpsiMass");

  string reducestr;

  TCanvas c1("c1","c1",10,10,1100,1000);
  c1.Divide(3,4);

  TCanvas c2("c2","c2",10,10,1100,1000);
  c2.Divide(3,4);
  
  string catstring_type;

  if(DataCat == 0) catstring_type = "GG";
  else if(DataCat == 1) catstring_type = "GT";
  else if(DataCat == 2) catstring_type = "TT";

  for(int i = 1;i<=nbinspt;i++){
    for(int j = 1;j<=nbinseta;j++){

       // ~ empty bins here
      if (DataCat == 0 && i == 1 && j < 3) continue;
      if (DataCat > 0 && i == 11) continue;

      RooPlot *mframe = JpsiMass->frame();

      ostringstream ostTitle;
      ostTitle << ptbinlimits[i-1] << " < pT < " << ptbinlimits[i] << " GeV, " << etabinlimits[j-1] << " < |#eta| < " << etabinlimits[j];

      reducestr = "Mass fit " + catstring_type + " muons: " + ostTitle.str();

      mframe->SetTitle(reducestr.c_str());

      string cutstring;
      string halfcutstring = nameCategoryCut(i,j);

      cutstring = halfcutstring + " && JpsiType == JpsiType::" + catstring_type;

      data->plotOn(mframe,DataError(RooAbsData::SumW2),MarkerSize(0.5),Cut(cutstring.c_str()));
    
      string catstring_pt = nameCategory("P",i);
      string catstring_eta = nameCategory("E",j);

      JpsiPtType->setLabel(catstring_pt.c_str());
      JpsiEtaType->setLabel(catstring_eta.c_str());

      string bkgname = nameSplitVar("expFunct",i,j);
      string totPDFname = nameSplitVar("totPDF",i,j);

      ws->pdf(totPDFname.c_str())->plotOn(mframe,Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));
      ws->pdf(totPDFname.c_str())->plotOn(mframe,DrawOption("F"),FillColor(kGreen),Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));
      ws->pdf(totPDFname.c_str())->plotOn(mframe,Components(*(ws->pdf(bkgname.c_str()))),DrawOption("F"),FillColor(kRed),Normalization(totalPDF->expectedEvents(RooArgSet(*JpsiMass)),RooAbsReal::NumEvent));

      data->plotOn(mframe,DataError(RooAbsData::SumW2),MarkerSize(0.5),Cut(cutstring.c_str()));

      if(j == 1) {c1.cd(i);mframe->Draw();}
      if(j == 2) {c2.cd(i);mframe->Draw();}

    }
  }

  reducestr = catstring_type + "massfit_barrel.gif";   
  c1.SaveAs(reducestr.c_str());

  reducestr = catstring_type + "massfit_endcap.gif";
  c2.SaveAs(reducestr.c_str());

  return;

}

void printResults(RooWorkspace *ws, const int DataCat)
{
  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooCategory *JpsiPtType = ws->cat("JpsiPtType");
  RooCategory *JpsiEtaType = ws->cat("JpsiEtaType");

  TCanvas c3("c3","c3",10,10,500,600);
  TCanvas c4("c4","c4",10,10,500,600);

  double ptbincenters[nbinspt] = {nbinspt*20.};
  double ptbinerrors[nbinspt] = {nbinspt*0.};
  double ycenters1[nbinspt] = {nbinspt*0.};
  double yerrors1[nbinspt] = {nbinspt*0.};
  double ycenters2[nbinspt] = {nbinspt*0.};
  double yerrors2[nbinspt] = {nbinspt*0.};

  string catstring_type;

  if(DataCat == 0) catstring_type = "GG";
  else if(DataCat == 1) catstring_type = "GT";
  else if(DataCat == 2) catstring_type = "TT";

  cout << endl << endl << catstring_type << " J/psi yields:" << endl;

  RooDataSet* reddata2; 
  string cutstring;

  for(int i = 1;i<=nbinspt;i++){

    ptbincenters[i-1] = (ptbinlimits[i] + ptbinlimits[i-1])/2. ;
    ptbinerrors[i-1] = ptbinlimits[i] - ptbincenters[i-1] ;

    for(int j = 1;j<=nbinseta;j++){

      // ~ empty bins here, if any
      if (DataCat == 0 && i == 1 && j < 3) continue;
      if (DataCat > 0 && i == 11) continue;

      string Nsigname = nameSplitVar("NSig",i,j);
      // string coeffGaussname = nameSplitVar("coeffGauss",i,j);
      // string sigmaSig1name = nameSplitVar("sigmaSig1",-1,j);
      // string sigmaSig2name = "sigmaSig2";

      string halfcutstring = nameCategoryCut(i,j);

      cutstring = "(MCType == MCType::PR || MCType == MCType::NP) && JpsiType == JpsiType::" + catstring_type + " && " + halfcutstring;

      const double Nsig   = ws->var(Nsigname.c_str())->getVal();
      const double errSig = ws->var(Nsigname.c_str())->getError();
     
      cout << "Fit for bin (" << i << "," << j << "): " << Nsig << " +/- " << errSig << endl;
      reddata2 = (RooDataSet*)data->reduce(RooArgList(*JpsiPtType),cutstring.c_str()); 
      cout << "True MC for bin (" << i << "," << j << "): " << reddata2->sumEntries() << endl;
      
      if (j == 1) {
	ycenters1[i-1] = Nsig - reddata2->sumEntries();
        yerrors1[i-1] = errSig;
      }
      if (j == 2) {
	ycenters2[i-1] = Nsig - reddata2->sumEntries();
        yerrors2[i-1] = errSig;
      }
    }
  }

  TGraphErrors *gpull1 = new TGraphErrors(nbinspt,ptbincenters,ycenters1,ptbinerrors,yerrors1);

  c3.cd();
  gpull1->SetTitle("Pull of fit results - barrel");
  gpull1->SetMarkerStyle(20);
  gpull1->SetMarkerColor(kRed);
  gpull1->Draw("AP");

  cutstring = catstring_type + "pull_barrel.gif";
  c3.SaveAs(cutstring.c_str());

  TGraphErrors *gpull2 = new TGraphErrors(nbinspt,ptbincenters,ycenters2,ptbinerrors,yerrors2);

  c4.cd();
  gpull2->SetTitle("Pull of fit results - endcap");
  gpull2->SetMarkerStyle(20);
  gpull2->SetMarkerColor(kBlue);
  gpull2->Draw("AP");

  cutstring = catstring_type + "pull_endcap.gif";
  c4.SaveAs(cutstring.c_str());

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

  // read back pT-eta values
 
  ifstream fpt;    fpt.open("ptranges.txt");
  if (!fpt) { cout << "Error opening file ptranges.txt" << endl; assert(0); }

  Double_t min = 0., max = 0.;  
  Int_t i = 0;  

  //Read in the file and store back to array
  while( !fpt.eof() ){
    fpt >> min >> max;
    ptbinlimits[i] = min;   i++;
  }
  ptbinlimits[i-1] = max;

  i = 0;  

  ifstream feta;
  feta.open("etaranges.txt");
  if (!feta){ cout << "Error opening file etaranges.txt" << endl; assert(0); }

  while(!feta.eof()){
    feta >> min >> max;
    etabinlimits[i] = min;  i++;    
  }
  etabinlimits[i-1] = max;

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *reddata = (RooDataSet*)fIn.Get("data");

  reddata->setWeightVar("MCweight");

  ws->import(*reddata);

  setRanges(ws);

  ws->var("JpsiMass")->setBins(60);

  setRanges(ws);

  //DEFINE SIGNAL AND BACKGROUND
  defineSignal(ws);
  defineBackground(ws);

  // Total PDF (signal CB+Gauss)
  ws->factory("SUM::totPDF(NSig[5000.,10.,10000000.]*sigCBGauss,NBkg[2000.,10.,10000000.]*expFunct)");

  //DO SIMULTANEOUS FIT
  RooSimWSTool wst(*ws);
  wst.build("totPDF_sim","totPDF",
	    //  SplitParam("NSig,NBkg,coefExp,coeffGauss,sigmaSig1","JpsiPtType,JpsiEtaType"));
	    SplitParam("NSig,NBkg,coefExp,coeffGauss,sigmaSig1","JpsiPtType,JpsiEtaType"));
	    // SplitParam("alpha,enne","JpsiEtaType")); 
  
  //Make subsamples to be used later

  RooDataSet *GGdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  RooDataSet *GTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  RooDataSet *TTdata = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::TT");
  RooDataHist *GGdataBin = new RooDataHist("GGdataBin","GGdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdata);
  RooDataHist *GTdataBin = new RooDataHist("GTdataBin","GTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GTdata);
  RooDataHist *TTdataBin = new RooDataHist("TTdataBin","TTdataBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*TTdata);

  //GG CASE 
 
  // fix some parameters  
  // ws->var("alpha")->setConstant(kTRUE);  
  // ws->var("enne")->setConstant(kTRUE);  

  if (sidebandPrefit) prefitSideband(ws,0); 

  if(prefitSignalMass){ 

    RooDataSet *GGdataTr = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG && (MCType == MCType::PR || MCType == MCType::NP)");

    // RooDataSet *GGdataTr1 = (RooDataSet*)GGdataTr->reduce("JpsiEtaType == JpsiEtaType::E1");
    // RooDataSet *GGdataTr2 = (RooDataSet*)GGdataTr->reduce("JpsiEtaType == JpsiEtaType::E2");
 
    RooDataHist *GGdataTrBin = new RooDataHist("GGdataTrBin","GGdataTrBin",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdataTr);

    // RooDataHist *GGdataTrBin1 = new RooDataHist("GGdataTrBin1","GGdataTrBin1",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdataTr1);

    // RooDataHist *GGdataTrBin2 = new RooDataHist("GGdataTrBin2","GGdataTrBin2",RooArgSet(*(ws->var("JpsiMass")),*(ws->cat("MCType")),*(ws->cat("JpsiPtType")),*(ws->cat("JpsiEtaType"))),*GGdataTr2);

    ws->pdf("sigCBGauss")->fitTo(*GGdataTrBin,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);

    /* // ws->pdf("sigCBGauss")->fitTo(*GGdataTr1,SumW2Error(kTRUE)); 
    ws->pdf("sigCBGauss")->fitTo(*GGdataTrBin1,SumW2Error(kTRUE));

    double enneVal = ws->var("enne")->getVal();
    ws->var("enne_E1")->setVal(enneVal);
    ws->var("enne_E1")->setConstant(kTRUE);

    double alphaVal = ws->var("alpha")->getVal();
    ws->var("alpha_E1")->setVal(alphaVal);
    ws->var("alpha_E1")->setConstant(kTRUE);

    // ws->pdf("sigCBGauss")->fitTo(*GGdataTr2,SumW2Error(kTRUE)); 
    ws->pdf("sigCBGauss")->fitTo(*GGdataTrBin2,SumW2Error(kTRUE));

    enneVal = ws->var("enne")->getVal();
    ws->var("enne_E2")->setVal(enneVal);
    ws->var("enne_E2")->setConstant(kTRUE);

    alphaVal = ws->var("alpha")->getVal();
    ws->var("alpha_E2")->setVal(alphaVal);
    ws->var("alpha_E2")->setConstant(kTRUE);*/
  
  }

  // ws->pdf("totPDF_sim")->fitTo(*GGdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF_sim")->fitTo(*GGdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,0);

  printResults(ws,0);

  //GT case

  // fix some parameters 
  /* ws->var("alpha_E1")->setConstant(kTRUE);
  ws->var("alpha_E2")->setConstant(kTRUE);
  ws->var("enne_E1")->setConstant(kTRUE); 
  ws->var("enne_E2")->setConstant(kTRUE); */

  if (sidebandPrefit) prefitSideband(ws,1);

  // ws->pdf("totPDF_sim")->fitTo(*GTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  ws->pdf("totPDF_sim")->fitTo(*GTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,1);

  printResults(ws,1);

  //TT case

  // fix some parameters 
  //ws->var("alpha")->setConstant(kTRUE); 
  //ws->var("enne")->setConstant(kTRUE); 

  /* if (sidebandPrefit) prefitSideband(ws,2);

  // ws->pdf("totPDF_sim")->fitTo(*TTdata,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));
  // ws->pdf("totPDF_sim")->fitTo(*TTdataBin,Extended(1),Save(1),Minos(0),NumCPU(2),SumW2Error(kTRUE));

  drawResults(ws,2);

  printResults(ws,2); */

  return 1;
}

