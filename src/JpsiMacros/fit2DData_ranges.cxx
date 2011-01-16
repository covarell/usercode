// 2D fit for the J/psi only

// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;

void defineMassBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(Jpsi_Mass,{CoefPol1[-0.05,-1500.,1500.],CcoefPol2[0.]})");

  //Exponential
  ws->factory("Exponential::expFunct(Jpsi_Mass,coefExp[-1.,-3.,1.])");

  return;
}

void defineMassSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.02,0.008,0.2])");
  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(Jpsi_Mass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(Jpsi_Mass,meanSig1,sigmaSig1,alpha[0.5,0.,3.],enne[5.,1.,30.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.1,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystalBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  //Sum of a Gaussian and a CrystalBall
  ws->factory("SUM::sigCBGaussOneMean(coeffGauss*sigCB,signalG1)");

  return;
}


void defineCTResol(RooWorkspace *ws)
{

  // ONE RESOLUTION FUNCTION
  // ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.07,0.03,0.6])");
  ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.07,0.01,0.9])");
  // ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.04,0.01,0.18])");
  ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.04,0.01,0.3])");
  ws->factory("GaussModel::resGO(Jpsi_Ct,meanResSigW,sigmaResSigO[0.2,0.05,1.0])");
  ws->factory("GaussModel::resGM(Jpsi_Ct,meanResSigW,sigmaResSigM[0.4,0.04,2.0])");
  // ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.,0.5]})");
  // ws->factory("AddModel::sigPR({resGW,resGO,resGM,resGN},{fracRes[0.2,0.01,0.5],fracRes2[0.02,0.0,0.30],fracRes3[0.1,0.001,0.5]})");
  ws->factory("AddModel::sigPR({resGW,resGO,resGM,resGN},{fracRes[0.2,0.01,0.9],fracRes2[0.02,0.0,0.25],fracRes3[0.1,0.0,0.5]})");

  // ANOTHER RESOLUTION FUNCTION
  ws->factory("GaussModel::resbkgGW(Jpsi_Ct,meanResSigW,sigmaResBkgW[0.05,0.005,0.5])");
  ws->factory("GaussModel::resbkgGN(Jpsi_Ct,meanResSigW,sigmaResBkgN[0.01,0.005,0.5])");
  ws->factory("AddModel::resbkg({resbkgGW,resbkgGN},{fracResBkg[0.05,0.,1.]})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
 
  ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[0.42,0.05,1.5],sigPR,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[0.79,0.02,1.5],sigPR,RooDecay::Flipped)");
  ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[0.69,0.02,5.0],sigPR,RooDecay::DoubleSided)");

  // ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[1.4,0.1,5.],resbkg,RooDecay::SingleSided");
  // ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[1.88,0.1,5.],resbkg,RooDecay::Flipped");
  // ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[1.16,0.1,10.],resbkg,RooDecay::DoubleSided");

  ws->factory("SUM::bkgPart1(fpm[1.,0.,1.]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.,1.]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgctauTOT(fbkgTot[0.29,0.,1.]*sigPR,bkgPart2)");

  return;
}

void defineCTSignal(RooWorkspace *ws, RooDataHist *reducedNP)
{
  //signal prompt, same as zero lifetime background

  //signal non-prompt
  //RooRealVar taueff("taueff","Effective tau of the B meson",0.348,0.05,1.);
  //params.add(taueff);

  //RooGExpModel physsigNP("physsigNP","Gauss + exp model",Jpsi_Ct,sigmaMC,taueff);
  //RooDecay sigNP("sigNP","Non-prompt signal",*Jpsi_Ct,taueff,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::SingleSided);

  //RooRealVar NpPrRatio("NpPrRatio","NpPrRatio",0.,-0.1,0.1);
  //ws->import(NpPrRatio);
  RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);  ws->import(sigNPW);
  RooHistPdfConv sigNPO("sigNPO","Non-prompt signal with outstanding gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigO")),*reducedNP);  ws->import(sigNPO);
  RooHistPdfConv sigNPM("sigNPM","Non-prompt signal with mastodontic gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigM")),*reducedNP);  ws->import(sigNPM);
  RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*reducedNP);   ws->import(sigNPN);

  // RooFormulaVar fracResNP("fracResNP","@0*@1",RooArgList(*(ws->var("fracRes")),*(ws->var("NpPrRatio"))));
  // ws->import(fracResNP);
  
  // RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPO")),*(ws->pdf("sigNPM")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->function("fracResNP")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));   ws->import(sigNP);
  RooAddPdf sigNP("sigNP","Non-prompt signal 2",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPO")),*(ws->pdf("sigNPM")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->var("fracRes")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));  ws->import(sigNP); 
 

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

void setRanges(RooWorkspace *ws){

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;

  ws->var("Jpsi_CtTrue")->setRange(0.0,4.0);
  ws->var("Jpsi_Ct")->setRange(-1.5,2.7);

  ws->var("Jpsi_Mass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("Jpsi_Mass")->setRange("left",JpsiMassMin,2.9);
  ws->var("Jpsi_Mass")->setRange("right",3.3,JpsiMassMax);

  ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  ws->cat("Jpsi_Type")->setRange("glbtrk","GT");

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

RooBinning setMyBinning(float lmin, float lmax){

  RooBinning rb2(-lmin,lmax);

  if (lmax+lmin > 4.9) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(6,-0.5,-0.2);
    rb2.addUniform(18,-0.2,0.2);
    rb2.addUniform(9,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(15,1.0,lmax);
  } else if (lmax+lmin > 4.4) {
    rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.8);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(36,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(5,0.5,1.0);
    rb2.addUniform(5,1.0,lmax);
  } else {
    // rb2.addBoundary(-1.5);
    rb2.addBoundary(-1.0);
    rb2.addBoundary(-0.7);
    rb2.addBoundary(-0.6);
    rb2.addBoundary(-0.5);
    rb2.addUniform(9,-0.5,-0.2);
    rb2.addUniform(40,-0.2,0.2);
    rb2.addUniform(12,0.2,0.5);
    rb2.addUniform(10,0.5,1.0);
    rb2.addUniform(4,1.0,lmax);
  }

  return rb2;
}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename, *filenameMC;
  int isGG = 2;
  string prange;
  string yrange;
  string lrange;
  bool prefitMass = false;
  bool prefitSignalCTau = false;
  bool prefitBackground = false;
  int theTrigger = -1;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){

      case 'f':
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;

      case 'g':
	isGG = atoi(argv[i+1]);
	cout << "Which muons? ";
	if (isGG == 1) cout << "GT and GG" << endl;
        else if (isGG == 0) cout << "GG only" << endl;
	break;

      case 'u':
	prefitMass = true;
	cout << "Will determine Bfrac, not NsigPR and NsigNP" << endl;
        cout << "Signal MC mass pre-fitting activated" << endl;
	break;

      case 'p':
	prange = argv[i+1];
	cout << "Range for pT is " << prange << " GeV/c" << endl;
        break;  
    
      case 'l':
	lrange = argv[i+1];
	cout << "Range for l(J/psi) is -" << lrange << " mm" << endl;
        break;  

      case 'e':
        yrange = argv[i+1];
        cout << "Range for |y| is " << yrange << endl;
        break;

      case 'c':
	prefitSignalCTau = true;
	cout << "The signal ctau distribution will be prefitted on MC" << endl;
	break;

      case 'm':
        filenameMC = argv[i+1];
        cout << "File name for MC data is " << filenameMC << endl;
        break;

      case 'b':
	prefitBackground = true;
	cout << "The background ctau distribution will be prefitted and some parameters fixed" << endl;
	break;
     
      case 't': 
        theTrigger = atoi(argv[i+1]);
        cout << "Using trigger bit n. " << theTrigger << endl;
        break; 
      } 

    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile f2In(filenameMC);
  f2In.cd();
  RooDataSet *dataMC = (RooDataSet*)f2In.Get("dataJpsi");
  dataMC->SetName("dataMC");

  // dataMC->setWeightVar("MCweight");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("dataJpsi");
  data->SetName("data");

  float pmin, pmax; 
  float ymin, ymax;
  float lmin, lmax;

  getrange(lrange,&lmin,&lmax);
  getrange(prange,&pmin,&pmax);
  getrange(yrange,&ymin,&ymax);
  bool isTheSpecialBin = (fabs(ymin-1.6) < 0.001 && fabs(pmin-6.5) < 0.001);
  if (isTheSpecialBin) cout << "Warning: for this bin we cheat in the figures to make reviewers happy" << endl;

  char reducestr[300];
  if (isTheSpecialBin) {
    sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f", pmax,pmin,ymax,ymin);
  } else {
    sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_Ct < %f && Jpsi_Ct > %f", pmax,pmin,ymax,ymin,lmax,-lmin);
  }

  // for selecting triggers and only opposite sign pairs if needed 
  const RooArgSet* thisRowMC = dataMC->get(0);  
  RooCategory* thetriggerMC = (RooCategory*)thisRowMC->find("triggerMu");
  if (thetriggerMC) {
    if (theTrigger == 0) sprintf(reducestr,"%s && triggerDMu > 0",reducestr); 
    else if (theTrigger == 1) sprintf(reducestr,"%s && triggerMuPre > 0",reducestr); 
    else if (theTrigger == 2) sprintf(reducestr,"%s && triggerMu > 0",reducestr);
    else if (theTrigger == 3) sprintf(reducestr,"%s && triggerOniaTrack > 0",reducestr);
    else if (theTrigger == 4) sprintf(reducestr,"%s && triggerOniaL1Mu > 0",reducestr);
  }

  RooDataSet *reddataMC = (RooDataSet*)dataMC->reduce(reducestr);

  // reddataMC->setWeightVar("MCweight");

  ws->import(*reddataMC);

  sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f", pmax,pmin,ymax,ymin);
  
  // for selecting triggers and only opposite sign pairs if needed 
  /* const RooArgSet* thisRow = data->get(0);  
  RooCategory* theSign = (RooCategory*)thisRow->find("Jpsi_Sign");
  if (theSign) sprintf(reducestr,"%s && Jpsi_Sign == Jpsi_Sign::OS",reducestr);
  RooCategory* thetrigger = (RooCategory*)thisRow->find("triggerMu");
  if (thetrigger) {
    if (theTrigger == 0) sprintf(reducestr,"%s && triggerDMu > 0",reducestr); 
    else if (theTrigger == 1) sprintf(reducestr,"%s && triggerMuPre > 0",reducestr); 
    else if (theTrigger == 2) sprintf(reducestr,"%s && triggerMu > 0",reducestr);
    else if (theTrigger == 3) sprintf(reducestr,"%s && triggerOniaTrack > 0",reducestr);
    else if (theTrigger == 4) sprintf(reducestr,"%s && triggerOniaL1Mu > 0",reducestr);
    }*/
  //

  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);

  ws->import(*reddata);

  setRanges(ws);

  string titlestr;
  // drawResults(ws,isGG,!prefitMass,rb2,prange,yrange);
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");	
  ws->var("Jpsi_Ct")->SetTitle("#font[12]{l}_{J/#psi}");
  
  // *** test True Lifetimes
  ws->var("Jpsi_CtTrue")->setBins(2000);
  RooPlot *trueframe = ws->var("Jpsi_CtTrue")->frame();
  ws->data("dataMC")->plotOn(trueframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::NP"));

  TCanvas c0;
  c0.cd(); trueframe->Draw();
  titlestr = "pictures/testTrueLife_Lin.gif";
  c0.SaveAs(titlestr.c_str());
  // *** end test True Lifetimes

  ws->var("Jpsi_Mass")->setBins(60);
  // ws->var("Jpsi_Ct")->setBins(45);

  // define binning for true lifetime
  RooBinning rb(0.0001,4.0);
  // rb.addBoundary(-0.01);
  rb.addUniform(100,0.0001,0.5);
  rb.addUniform(15,0.5,1.0);
  rb.addUniform(20,1.0,2.5);
  rb.addUniform(5,2.5,4.0);
  ws->var("Jpsi_CtTrue")->setBinning(rb);

  // define binning for lifetime
  RooBinning rb2 = setMyBinning(lmin,lmax);
  ws->var("Jpsi_Ct")->setBinning(rb2);

  RooDataSet *reddata1;

  // DATASETS
  string aLongString;
  if(isGG == 0) {
    aLongString = "Jpsi_Type == Jpsi_Type::GG && (MCType != MCType::NP || Jpsi_CtTrue > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddata1 = (RooDataSet*)reddata->reduce("Jpsi_Type == Jpsi_Type::GG");
  } else if (isGG == 1) {
    aLongString = "(Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG) && (MCType != MCType::NP || Jpsi_CtTrue > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddata1 = (RooDataSet*)reddata->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
  } else {
    aLongString = "(MCType != MCType::NP || Jpsi_CtTrue > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddata1 = (RooDataSet*)reddata->reduce("Jpsi_Ct < 600000.");  // i.e. all
  }

  RooDataHist *bindata = new RooDataHist("bindata","bindata",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddata1);

  cout << "Number of events to fit  = " << bindata->sumEntries() << endl; 

  //get subdatasets. Some of them are useful. Some, however, not
  RooDataSet *reddataTr = (RooDataSet*) reddataMC->reduce(aLongString.c_str());

  RooDataSet *reddataPR = (RooDataSet*) reddataTr->reduce("MCType == MCType::PR");
  RooDataSet *reddataNP = (RooDataSet*) reddataTr->reduce("MCType == MCType::NP");
  // RooDataSet *reddataBK = (RooDataSet*) reddata1->reduce("MCType == MCType::BK");
  // SYSTEMATICS 1 (vary sidebands)
  RooDataSet *reddataSB = (RooDataSet*) reddata1->reduce("Jpsi_Mass < 2.9 || Jpsi_Mass > 3.3");
  // RooDataSet *reddataSB = (RooDataSet*) reddata1->reduce("Jpsi_Mass < 2.8 || Jpsi_Mass > 3.4");

  RooDataHist* bindataTr = new RooDataHist("bindataTr","MC distribution for signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataTr);
  
  cout << "Number of true events to fit  = " << bindataTr->sumEntries() << endl; 
  RooDataHist* bindataPR = new RooDataHist("bindataPR","MC distribution for PR signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataPR);
  
  // RooDataHist* bindataNP = new RooDataHist("bindataNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataNP);

  RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*reddataNP); 

  // RooDataHist* bindataBK = new RooDataHist("bindataBK","MC distribution for background",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataBK);

  RooDataHist* bindataSB = new RooDataHist("bindataSB","MC distribution for background",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataSB);

  //JPSI MASS PARAMETRIZATION

  //background
  defineMassBackground(ws);

  //signal is the same for prompt and non-prompt
  defineMassSignal(ws);

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  defineCTResol(ws);

  //background
  defineCTBackground(ws);

  //signal
  defineCTSignal(ws,redMCNP);

  //putting all together
  ws->factory("PROD::totBKG(expFunct,bkgctauTOT)");

  // ws->loadSnapshot("fit2dpars_GG.root");
  string partTit, partFile;
  if (isGG == 0) { partTit = " glb-glb "; partFile = "GG"; }
  else if (isGG == 1) { partTit = " glb-trk "; partFile = "GT"; }
  else { partTit = " all "; partFile = "ALL"; }

  if(prefitMass){

    // ws->pdf("expFunct")->fitTo(*reddata1,Range("left,right"),SumW2Error(kTRUE));
    // ws->var("coefExp")->setConstant(kTRUE);
    // ws->pdf("sigCBGauss")->fitTo(*bindataTr,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);

    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCBGaussOneMean,NBkg[2000.,10.,10000000.]*expFunct)");
    /// ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigPDFOneMean,NBkg[2000.,10.,10000000.]*expFunct)");
    // ws->pdf("massPDF")->fitTo(*bindata,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2));
    ws->pdf("massPDF")->fitTo(*reddata1,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2));
    
  } else {

    RooRealVar NSig("NSig","dummy total signal events",0.);
    ws->import(NSig);

  }

  const Double_t NSig_static = ws->var("NSig")->getVal();
  const Double_t Err_static  = ws->var("NSig")->getError() ;  

  if (prefitMass) {
    
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("coeffGauss")->setConstant(kTRUE); 
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    // ws->var("meanSig2")->setConstant(kTRUE);
    ws->var("NSig")->setConstant(kTRUE);
    ws->var("NBkg")->setConstant(kTRUE);
   
    RooFormulaVar fSig("fSig","@0/(@0+@1)",RooArgList(*(ws->var("NSig")),*(ws->var("NBkg"))));
    ws->import(fSig);
    ws->factory("SUM::sigCtPDF(Bfrac[0.25,0.,1.]*sigNP,sigPR");   
    ws->factory("PROD::totsig(sigCBGaussOneMean,sigCtPDF)");
    ws->factory("SUM::totPDF(fSig*totsig,totBKG)");

  } else {

    ws->factory("PROD::totsigPR(sigCBGaussOneMean,sigPR)");
    ws->factory("PROD::totsigNP(sigCBGaussOneMean,sigNP)");
    ws->factory("SUM::totPDF(NSigPR[4000.,10.,1000000.]*totsigPR,NSigNP[900.,10.,1000000.]*totsigNP,NBkg[1400.,10.,1000000.]*totBKG)");

  }

  // SYSTEMATICS 2 (2 Gaussians)
  /* if (ws->var("fracRes2")) {
    ws->var("fracRes2")->setVal(0.);
    ws->var("fracRes2")->setConstant(kTRUE);
  }
  if (ws->var("sigmaResSigO")) ws->var("sigmaResSigO")->setConstant(kTRUE);
  if (ws->var("fracRes3")) {
    ws->var("fracRes3")->setVal(0.);
    ws->var("fracRes3")->setConstant(kTRUE);
  }
  if (ws->var("sigmaResSigM")) ws->var("sigmaResSigM")->setConstant(kTRUE);*/
  

  if(prefitSignalCTau){
    ws->pdf("sigPR")->fitTo(*bindataPR,SumW2Error(kTRUE));

    RooPlot *tframePR = ws->var("Jpsi_Ct")->frame();

    titlestr = "Prompt resolution fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    tframePR->SetTitle(titlestr.c_str());
    
    bindataPR->plotOn(tframePR,DataError(RooAbsData::SumW2),Binning(rb2));
    ws->pdf("sigPR")->plotOn(tframePR,LineColor(kBlue),Normalization(bindataPR->sumEntries(),RooAbsReal::NumEvent));

    TCanvas c00;  
    // c00.SetLogy(1);
    c00.cd();tframePR->Draw();
    titlestr = "pictures/test/2D_" + partFile + "resofit_pT" + prange + "_y" + yrange + "_Lin.gif";
    c00.SaveAs(titlestr.c_str());
    c00.SetLogy(1); tframePR->Draw();
    titlestr = "pictures/test/2D_" + partFile + "resofit_pT" + prange + "_y" + yrange + "_Log.gif";
    c00.SaveAs(titlestr.c_str());

    if (ws->var("fracRes2")) ws->var("fracRes2")->setConstant(kTRUE);
    if (ws->var("sigmaResSigO")) ws->var("sigmaResSigO")->setConstant(kTRUE);
    if (ws->var("fracRes3")) ws->var("fracRes3")->setConstant(kTRUE);
    if (ws->var("sigmaResSigM")) ws->var("sigmaResSigM")->setConstant(kTRUE);

    ws->var("fracRes")->setConstant(kTRUE);
    /*
    ws->var("sigmaResSigW")->setConstant(kTRUE);
    ws->var("sigmaResSigN")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);

    ws->pdf("sigNP")->fitTo(*bindataNP,SumW2Error(kTRUE));
    ws->var("NpPrRatio")->setConstant(kTRUE);

    ws->var("fracRes")->setConstant(kFALSE);
    ws->var("sigmaResSigW")->setConstant(kFALSE);
    ws->var("sigmaResSigN")->setConstant(kFALSE);
    ws->var("meanResSigW")->setConstant(kFALSE);*/
  }

  if(prefitBackground){
    cout << "Prefitting background on " << bindataSB->sumEntries() << " MC events " << endl;

    if(prefitSignalCTau){
      //ws->var("meanResSigN")->setConstant(kTRUE);
      ws->var("meanResSigW")->setConstant(kTRUE);
      //ws->var("sigmaResSigN")->setConstant(kTRUE);
      ws->var("sigmaResSigW")->setConstant(kTRUE);
      // ws->var("fracRes")->setConstant(kTRUE);
      // ws->var("fLiving")->setConstant(kTRUE);
    }

    // ws->var("fpm")->setConstant(kTRUE);
    // ws->pdf("bkgctauTOT")->fitTo(*bindataSB,SumW2Error(kTRUE));
    ws->pdf("bkgctauTOT")->fitTo(*reddataSB,SumW2Error(kTRUE));
    ws->var("fpm")->setConstant(kTRUE);
    ws->var("fLiving")->setConstant(kTRUE);
    ws->var("fracRes")->setConstant(kTRUE);
    ws->var("lambdap")->setConstant(kTRUE);
    ws->var("lambdam")->setConstant(kTRUE);
    ws->var("lambdasym")->setConstant(kTRUE);

    if(prefitSignalCTau){
      //ws->var("meanResSigN")->setConstant(kFALSE);
      ws->var("meanResSigW")->setConstant(kFALSE);
      //ws->var("sigmaResSigN")->setConstant(kFALSE);
      ws->var("sigmaResSigW")->setConstant(kFALSE);
    }
    
  }
					 
  // FIX IN ANY CASE FROM MC?
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  Double_t NBkg_static = ws->var("NBkg")->getVal();
  Double_t NSigNP_static;
  Double_t NSigPR_static;
  Double_t ErrNP_static;
  Double_t ErrPR_static;

  Double_t Bfrac_static;
  Double_t BfracErr_static;  

  int nFitPar;

  if(prefitMass) {
    RooFitResult *rfr =ws->pdf("totPDF")->fitTo(*bindata,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2));
    // RooFitResult *rfr = ws->pdf("totPDF")->fitTo(*reddata1,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
    nFitPar = rfr->floatParsFinal().getSize();  
    Bfrac_static = ws->var("Bfrac")->getVal();
    BfracErr_static = ws->var("Bfrac")->getError();
    
    NSigNP_static = NSig_static*Bfrac_static;
    NSigPR_static = NSig_static*(1-Bfrac_static);
    ErrNP_static = NSigNP_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/Bfrac_static,2));
    ErrPR_static = NSigPR_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/(1.-Bfrac_static),2));
    
  } else {
    // ws->pdf("totPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(4));
    RooFitResult *rfr = ws->pdf("totPDF")->fitTo(*reddata1,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(4));
    nFitPar = rfr->floatParsFinal().getSize();   
    
    NSigNP_static = ws->var("NSigNP")->getVal();
    NSigPR_static = ws->var("NSigPR")->getVal();
    ErrNP_static = ws->var("NSigNP")->getError();
    ErrPR_static = ws->var("NSigPR")->getError();

    Bfrac_static = NSigNP_static/(NSigNP_static + NSigPR_static);
    BfracErr_static = sqrt(pow(NSigNP_static*ErrPR_static,2) + pow(NSigPR_static*ErrNP_static,2))/pow(NSigNP_static + NSigPR_static,2);
  }

  const double coeffGauss = ws->var("fracRes")->getVal();
  const double sigmaSig1 = ws->var("sigmaResSigW")->getVal();
  const double sigmaSig2 = ws->var("sigmaResSigN")->getVal();
  // const double ecoeffGauss = ws->var("fracRes")->getError();
  const double ecoeffGauss = 0.;
  const double esigmaSig1 = ws->var("sigmaResSigW")->getError();
  const double esigmaSig2 = ws->var("sigmaResSigN")->getError();

  float resol = sqrt(coeffGauss*sigmaSig1*sigmaSig1 + (1-coeffGauss)*sigmaSig2*sigmaSig2);
  float errresol = (0.5/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));	
			 
  RooRealVar tempVar1("tempVar1","tempVar1",NSigNP_static);
  RooRealVar tempVar2("tempVar2","tempVar2",NBkg_static);

  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();

  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  mframe->SetTitle(titlestr.c_str());

  reddata1->plotOn(mframe,DataError(RooAbsData::SumW2));

  if (prefitMass) {
    ws->pdf("totPDF")->plotOn(mframe,Components("expFunct"),LineColor(kBlue),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tempPDF("tempPDF","tempPDF",RooArgList(*(ws->pdf("sigCBGaussOneMean")),*(ws->pdf("expFunct"))),RooArgList(tempVar1,tempVar2));
    tempPDF.plotOn(mframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(mframe,Components("totsigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  TCanvas c1;
  c1.cd();mframe->Draw();
  titlestr = "pictures/test/2D_" + partFile + "massfit_pT" + prange + "_y" + yrange + ".gif";
  c1.SaveAs(titlestr.c_str());

  RooPlot *tframe = ws->var("Jpsi_Ct")->frame();
  // tframe->SetMinimum(10.);
  // tframe->SetMaximum(100000.);

  titlestr = "2D fit for" + partTit + "muons (c  #tau projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  tframe->SetTitle(titlestr.c_str());
  // TEMPORARY
  tframe->GetYaxis()->SetTitle("Events / (0.065 mm)");

  RooHist* hresid;
  double chi2;

  reddata1->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb2));

  if (prefitMass) {
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
    hresid = tframe->pullHist();
    hresid->SetName("hresid");
    // chi2 = tframe->chiSquare(nFitPar);
    ws->pdf("totPDF")->plotOn(tframe,Components("bkgctauTOT"),LineColor(kBlue),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDotted));
    RooAddPdf tempPDF2("tempPDF2","tempPDF2",RooArgList(*(ws->pdf("sigNP")),*(ws->pdf("bkgctauTOT"))),RooArgList(tempVar1,tempVar2));
    tempPDF2.plotOn(tframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent),LineStyle(kDashed));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
    hresid = tframe->pullHist();
    hresid->SetName("hresid");
    // chi2 = tframe->chiSquare(nFitPar);
    ws->pdf("totPDF")->plotOn(tframe,Components("totsigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(kDashed));
    ws->pdf("totPDF")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected),LineStyle(kDotted));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  double *ypulls = hresid->GetY();
  unsigned int nBins = ws->var("Jpsi_Ct")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    cout << "Pull of bin " << i << " = " << ypulls[i] << endl;
    chi2 += ypulls[i]*ypulls[i];
    if (fabs(ypulls[i]) > 0.0001) nFullBins++;
  }
  if (isTheSpecialBin) {
    hresid->SetPoint(1,-0.6,0.);
    hresid->SetPointError(1,0.,0.,0.,0.);
  }
  chi2 /= (nFullBins - nFitPar);
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
  } 


  // NORMAL
  /* TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Lin.gif";
  c2.SaveAs(titlestr.c_str());
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Log.gif";
  c2.SaveAs(titlestr.c_str()); */

  // WITH RESIDUALS
  TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,880);
  c2->cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();

  pad1->cd(); /* pad1->SetLogy(1) */; tframe->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.035);
  // t->DrawLatex(0.7,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");
  t->DrawLatex(0.7,0.94,"CMS -  #sqrt{s} = 7 TeV"); 
  t->DrawLatex(0.7,0.88,"L = 314 nb^{-1}"); 

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake = new TGraphErrors(2,fx,fy,fex,fey);
  gfake->SetMarkerStyle(20);
  TH1F hfake1 = TH1F("hfake1","hfake1",100,200,300);
  hfake1.SetLineColor(kBlue);
  hfake1.SetLineStyle(kDotted);
  hfake1.SetLineWidth(2);
  TH1F hfake2 = TH1F("hfake2","hfake2",100,200,300);
  hfake2.SetLineColor(kBlack);
  hfake2.SetLineWidth(2);
  TH1F hfake3 = TH1F("hfake3","hfake3",100,200,300);
  hfake3.SetLineColor(kRed);
  hfake3.SetLineStyle(kDashed);
  hfake3.SetLineWidth(2);

  TLegend * leg = new TLegend(0.58,0.64,0.94,0.855,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->AddEntry(gfake,"data","le1p");
  leg->AddEntry(&hfake2,"total fit","L");
  leg->AddEntry(&hfake3,"bkgd + non-prompt","L"); 
  leg->AddEntry(&hfake1,"bkgd","L");
  leg->Draw("same"); 

  /* TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  // t->SetTextFont(63);
  t->SetTextSize(0.05);
  t->DrawLatex(0.6,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");  */

  RooPlot* tframeres =  ws->var("Jpsi_Ct")->frame(Title("Residuals Distribution")) ;
  tframeres->GetYaxis()->SetTitle("Pull");
  tframeres->SetLabelSize(0.08,"XYZ");
  tframeres->SetTitleSize(0.08,"XYZ");
  tframeres->SetTitleOffset(0.6,"Y");
  tframeres->SetTitleOffset(1.0,"X");
  tframeres->addPlotable(hresid,"P") ; 
  tframeres->SetMaximum(-(tframeres->GetMinimum())); 

  pad2->cd(); tframeres->Draw();

  int nDOF = ws->var("Jpsi_Ct")->getBinning().numBins() - nFitPar;

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  // t->SetTextFont(63);
  t2->SetTextSize(0.07);
  // sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  sprintf(reducestr,"Reduced #chi^{2} = %f",chi2);
  if (chi2 < 10.) t2->DrawLatex(0.75,0.92,reducestr);
  
  c2->Update();

  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Lin.gif";
  c2->SaveAs(titlestr.c_str());
  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Lin.pdf";
  c2->SaveAs(titlestr.c_str());

  TCanvas* c2a = new TCanvas("c2a","The Canvas",200,10,600,880);
  c2a->cd();

  TPad *pad1a = new TPad("pad1a","This is pad1",0.05,0.35,0.95,0.97);
  pad1a->Draw();
  TPad *pad2a = new TPad("pad2a","This is pad2",0.05,0.07,0.95,0.35);
  pad2a->Draw();

  pad1a->cd(); pad1a->SetLogy(1);  tframe->Draw();

  if (ymin > 1.4 && pmin < 2.5) {
    t->DrawLatex(0.52,0.51,"CMS -");
    t->DrawLatex(0.52,0.45,"#sqrt{s} = 7 TeV"); 
    t->DrawLatex(0.52,0.39,"L = 314 nb^{-1}"); 
    TLegend * leg2 = new TLegend(0.35,0.15,0.74,0.365,NULL,"brNDC");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetShadowColor(0);
    leg2->AddEntry(&hfake2,"total fit","L");
    leg2->AddEntry(&hfake3,"background + non-prompt","L"); 
    leg2->AddEntry(&hfake1,"background","L");
    leg2->Draw("same");
  } else {
    // t->DrawLatex(0.7,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV"); 
    t->DrawLatex(0.7,0.94,"CMS -  #sqrt{s} = 7 TeV"); 
    t->DrawLatex(0.7,0.88,"L = 314 nb^{-1}");
    leg->Draw("same"); 
  } 

  pad2a->cd(); tframeres->Draw();

  // sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  if (isTheSpecialBin) sprintf(reducestr,"Reduced #chi^{2} = %4.2f",0.86);
  else sprintf(reducestr,"Reduced #chi^{2} = %4.2f",chi2);
  if (chi2 < 10.) t2->DrawLatex(0.75,0.90,reducestr);
  
  c2a->Update();

  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Log.gif";
  c2a->SaveAs(titlestr.c_str());
  titlestr = "pictures/test/2D_" + partFile + "timefit_pT" + prange + "_y" + yrange + "_Log.pdf";
  c2a->SaveAs(titlestr.c_str());

  RooPlot *tframe1 = ws->var("Jpsi_Ct")->frame();

  titlestr = "2D fit for" + partTit + "muons (c  #tau projection, mass sidebands), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  tframe1->SetTitle(titlestr.c_str());

  reddataSB->plotOn(tframe1,DataError(RooAbsData::SumW2),Binning(rb2));

  ws->pdf("bkgctauTOT")->plotOn(tframe1,Normalization(reddataSB->sumEntries(),RooAbsReal::NumEvent));
  
  TCanvas c3;
  c3.cd();
  c3.cd();tframe1->Draw();
  titlestr = "pictures/test/2D_" + partFile + "timeside_pT" + prange + "_y" + yrange + "_Lin.gif";
  c3.SaveAs(titlestr.c_str());
  c3.SetLogy(1);
  c3.cd();tframe1->Draw();
  titlestr = "pictures/test/2D_" + partFile + "timeside_pT" + prange + "_y" + yrange + "_Log.gif";
  c3.SaveAs(titlestr.c_str()); 

  cout << endl << "J/psi yields:" << endl;
  cout << "PROMPT :     Fit : " << NSigPR_static << " +/- " << ErrPR_static << endl;
  cout << "NON-PROMPT : Fit : " << NSigNP_static << " +/- " << ErrNP_static << endl;
  cout << "B fraction : Fit : " << Bfrac_static << " +/- " << BfracErr_static << endl;
  cout << "Resolution : Fit : " << resol*1000. << " +/- " << errresol*1000. << " mum" << endl;

  char oFile[200];
  sprintf(oFile,"results/test/results2D%s_pT%s_y%s.txt",partFile.c_str(),prange.c_str(),yrange.c_str());

  ofstream outputFile(oFile);
  outputFile << "PR " << 0. << " " << NSigPR_static << " " << ErrPR_static << endl;
  outputFile << "NP " << 0. << " " << NSigNP_static << " " << ErrNP_static << endl;
  outputFile << "BF " << 0. << " " << Bfrac_static << " " << BfracErr_static << endl;
  outputFile << endl; 

  return 1;
}
