// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
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
  ws->factory("Polynomial::CPolFunct(JpsiMass,{CoefPol1[-0.05,-1500.,1500.],CcoefPol2[0.]})");

  //Exponential
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-1.,-3.,0.1])");

  return;
}

void defineMassSignal(RooWorkspace *ws)
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


void defineCTResol(RooWorkspace *ws)
{

  // ONE RESOLUTION FUNCTION
  ws->factory("GaussModel::resGW(Jpsict,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.02,0.008,0.2])");
  // ws->factory("GaussModel::resGN(Jpsict,meanResSigN[0.,-0.1,0.1],sigmaResSigN[0.01,0.005,0.5])");
  ws->factory("GaussModel::resGN(Jpsict,meanResSigW,sigmaResSigN[0.02,0.01,0.03])");
  ws->factory("GaussModel::resGO(Jpsict,meanResSigW,sigmaResSigO[0.2,0.1,0.3])");
  ws->factory("GaussModel::resGM(Jpsict,meanResSigW,sigmaResSigM[0.4,0.3,0.5])");
  // ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.,0.5]})");
  ws->factory("AddModel::sigPR({resGW,resGO,resGM,resGN},{fracRes[0.4,0.3,0.8],fracRes2[0.02,0.0,0.15],fracRes3[0.1,0.0,0.5]})");

  // ANOTHER RESOLUTION FUNCTION
  ws->factory("GaussModel::resbkgGW(Jpsict,meanResSigW,sigmaResBkgW[0.05,0.005,0.5])");
  // ws->factory("GaussModel::resGN(Jpsict,meanResSigN[0.,-0.1,0.1],sigmaResSigN[0.01,0.005,0.5])");
  ws->factory("GaussModel::resbkgGN(Jpsict,meanResSigW,sigmaResBkgN[0.01,0.005,0.5])");
  ws->factory("AddModel::resbkg({resbkgGW,resbkgGN},{fracResBkg[0.05,0.,1.]})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
 
  ws->factory("Decay::bkg2(Jpsict,lambdap[0.42,0.0001,2.],sigPR,RooDecay::SingleSided");
  ws->factory("Decay::bkg3(Jpsict,lambdam[0.79,0.0001,2.],sigPR,RooDecay::Flipped");
  ws->factory("Decay::bkg4(Jpsict,lambdasym[0.69,0.0001,2.],sigPR,RooDecay::DoubleSided");
  // ws->factory("Decay::bkg2(Jpsict,lambdap[1.4,0.1,5.],resbkg,RooDecay::SingleSided");
  // ws->factory("Decay::bkg3(Jpsict,lambdam[1.88,0.1,5.],resbkg,RooDecay::Flipped");
  // ws->factory("Decay::bkg4(Jpsict,lambdasym[1.16,0.1,10.],resbkg,RooDecay::DoubleSided");

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

  //RooGExpModel physsigNP("physsigNP","Gauss + exp model",Jpsict,sigmaMC,taueff);
  //RooDecay sigNP("sigNP","Non-prompt signal",*Jpsict,taueff,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::SingleSided);

  //RooRealVar NpPrRatio("NpPrRatio","NpPrRatio",0.,-0.1,0.1);
  //ws->import(NpPrRatio);
  RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);  ws->import(sigNPW);
  RooHistPdfConv sigNPO("sigNPO","Non-prompt signal with outstanding gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigO")),*reducedNP);  ws->import(sigNPO);
  RooHistPdfConv sigNPM("sigNPM","Non-prompt signal with mastodontic gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigM")),*reducedNP);  ws->import(sigNPM);
  RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*reducedNP);   ws->import(sigNPN);

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

  ws->var("JpsictTrue")->setRange(0.0,4.0);
  ws->var("Jpsict")->setRange(-1.0,3.5);

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

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename, *filenameMC;
  int isGG = 2;
  string prange;
  string etarange;
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
	cout << "Is global global? ";
	if (isGG == 1) cout << "No" << endl;
        else if (isGG == 0) cout << "Yes" << endl;
        else cout << "Maybe" << endl;
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
       
      case 'e':
        etarange = argv[i+1];
        cout << "Range for |eta| is " << etarange << endl;
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
  RooDataSet *dataMC = (RooDataSet*)f2In.Get("data");
  dataMC->SetName("dataMC");

  dataMC->setWeightVar("MCweight");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");

  float pmin, pmax; 
  float etamin, etamax;

  getrange(prange,&pmin,&pmax);
  getrange(etarange,&etamin,&etamax);

  char reducestr[300], reducestr2[300];
  sprintf(reducestr,"JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f", pmax,pmin,etamax,etamin);

  if (theTrigger == 0) sprintf(reducestr2,"%s && trigger0 > 0",reducestr);
  else if (theTrigger == 1) sprintf(reducestr2,"%s && trigger1 > 0",reducestr);

  RooDataSet *reddataMC = (RooDataSet*)dataMC->reduce(reducestr);

  reddataMC->setWeightVar("MCweight");

  ws->import(*reddataMC);
  
  // for selecting only opposite sign pairs if needed 
  const RooArgSet* thisRow = data->get(0);  
  RooCategory* theSign = (RooCategory*)thisRow->find("JpsiSign");
  if (theSign) sprintf(reducestr2,"%s && JpsiSign == JpsiSign::OS",reducestr); 
  //

  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);

  ws->import(*reddata);

  setRanges(ws);

  string titlestr;
  
  // *** test True Lifetimes
  ws->var("JpsictTrue")->setBins(2000);
  RooPlot *trueframe = ws->var("JpsictTrue")->frame();
  ws->data("dataMC")->plotOn(trueframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::NP"));

  TCanvas c0;
  c0.cd(); trueframe->Draw();
  titlestr = "pictures/testTrueLife_Lin.gif";
  c0.SaveAs(titlestr.c_str());
  // *** end test True Lifetimes

  ws->var("JpsiMass")->setBins(60);
  // ws->var("Jpsict")->setBins(45);

  // define binning
  RooBinning rb(0.0001,4.0);
  // rb.addBoundary(-0.01);
  rb.addUniform(100,0.0001,0.5);
  rb.addUniform(15,0.5,1.0);
  rb.addUniform(20,1.0,2.5);
  rb.addUniform(5,2.5,4.0);
  ws->var("JpsictTrue")->setBinning(rb);

  // define binning
  RooBinning rb2(-1.0,3.5);
  rb2.addBoundary(-0.5);
  // rb2.addBoundary(-0.2);
  // rb2.addBoundary(-0.1);
  // rb2.addBoundary(-0.01);
  rb2.addUniform(10,-0.5,-0.2);
  rb2.addUniform(20,-0.2,0.2);
  rb2.addUniform(10,0.2,0.5);
  rb2.addUniform(5,0.5,1.0);
  rb2.addUniform(10,1.0,3.5);
  ws->var("Jpsict")->setBinning(rb2);

  RooDataSet *reddata1, *reddataMC1;

  // DATASETS
  if(isGG == 0) {
    reddataMC1 = (RooDataSet*)reddataMC->reduce("JpsiType == JpsiType::GG && (MCType != MCType::NP || JpsictTrue > 0.0001)");
    reddata1 = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GG");
  } else if (isGG == 1) {
    reddataMC1 = (RooDataSet*)reddataMC->reduce("JpsiType == JpsiType::GT && (MCType != MCType::NP || JpsictTrue > 0.0001)");
    reddata1 = (RooDataSet*)reddata->reduce("JpsiType == JpsiType::GT");
  } else {
    reddataMC1 = (RooDataSet*)reddataMC->reduce("MCType != MCType::NP || JpsictTrue > 0.0001");
    reddata1 = (RooDataSet*)reddata->reduce("Jpsict < 600000.");  // i.e. all
  }
  
  RooDataHist *bindata = new RooDataHist("bindata","bindata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddata1);

  cout << "Number of events to fit  = " << bindata->sumEntries() << endl; 

  //get subdatasets. Some of them are useful. Some, however, not
  RooDataSet *reddataTr = (RooDataSet*) reddataMC1->reduce("MCType == MCType::PR || MCType == MCType::NP");
  RooDataSet *reddataPR = (RooDataSet*) reddataMC1->reduce("MCType == MCType::PR");
  RooDataSet *reddataNP = (RooDataSet*) reddataMC1->reduce("MCType == MCType::NP");
  // RooDataSet *reddataBK = (RooDataSet*) reddata1->reduce("MCType == MCType::BK");
  RooDataSet *reddataSB = (RooDataSet*) reddata1->reduce("JpsiMass < 2.9 || JpsiMass > 3.3");

  RooDataHist* bindataTr = new RooDataHist("bindataTr","MC distribution for signal",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataTr);
  
  cout << "Number of true events to fit  = " << bindataTr->sumEntries() << endl; 
  RooDataHist* bindataPR = new RooDataHist("bindataPR","MC distribution for PR signal",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataPR);
  
  RooDataHist* bindataNP = new RooDataHist("bindataNP","MC distribution for NP signal",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataNP);

  RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("JpsictTrue"))),*reddataNP); 

  // RooDataHist* bindataBK = new RooDataHist("bindataBK","MC distribution for background",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataBK);

  RooDataHist* bindataSB = new RooDataHist("bindataSB","MC distribution for background",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataSB);

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

  if(prefitMass){

    ws->pdf("sigCBGauss")->fitTo(*bindataTr,SumW2Error(kTRUE));
    ws->var("enne")->setConstant(kTRUE);

    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCBGauss,NBkg[2000.,10.,10000000.]*expFunct)");
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
    ws->var("coeffGauss")->setConstant(kTRUE); 
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("NSig")->setConstant(kTRUE);
    ws->var("NBkg")->setConstant(kTRUE);
   
    RooFormulaVar fSig("fSig","@0/(@0+@1)",RooArgList(*(ws->var("NSig")),*(ws->var("NBkg"))));
    ws->import(fSig);
    ws->factory("SUM::sigCtPDF(Bfrac[0.25,0.,1.]*sigNP,sigPR");   
    ws->factory("PROD::totsig(sigCBGauss,sigCtPDF)");
    ws->factory("SUM::totPDF(fSig*totsig,totBKG)");

  } else {

    ws->factory("PROD::totsigPR(sigCBGauss,sigPR)");
    ws->factory("PROD::totsigNP(sigCBGauss,sigNP)");
    ws->factory("SUM::totPDF(NSigPR[4000.,10.,1000000.]*totsigPR,NSigNP[900.,10.,1000000.]*totsigNP,NBkg[1400.,10.,1000000.]*totBKG)");

  }

  if(prefitSignalCTau){
    ws->pdf("sigPR")->fitTo(*bindataPR,SumW2Error(kTRUE));

    if (ws->var("fracRes2")) ws->var("fracRes2")->setConstant(kTRUE);
    if (ws->var("sigmaResSigO")) ws->var("sigmaResSigO")->setConstant(kTRUE);
    if (ws->var("fracRes3")) ws->var("fracRes3")->setConstant(kTRUE);
    if (ws->var("sigmaResSigM")) ws->var("sigmaResSigM")->setConstant(kTRUE);

    /* ws->var("fracRes")->setConstant(kTRUE);
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
    }

    ws->var("fpm")->setConstant(kTRUE);
    // ws->pdf("bkgctauTOT")->fitTo(*bindataSB,SumW2Error(kTRUE));
    ws->pdf("bkgctauTOT")->fitTo(*reddataSB,SumW2Error(kTRUE));
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

  if(prefitMass) {
    // ws->pdf("totPDF")->fitTo(*bindata,Minos(0),SumW2Error(kTRUE),NumCPU(4));
    ws->pdf("totPDF")->fitTo(*reddata,Minos(0),SumW2Error(kTRUE),NumCPU(4));

    Bfrac_static = ws->var("Bfrac")->getVal();
    BfracErr_static = ws->var("Bfrac")->getError();
    
    NSigNP_static = NSig_static*Bfrac_static;
    NSigPR_static = NSig_static*(1-Bfrac_static);
    ErrNP_static = NSigNP_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/Bfrac_static,2));
    ErrPR_static = NSigPR_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/(1.-Bfrac_static),2));
    
  } else {
    // ws->pdf("totPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(4));
    ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(4));

    // ws->saveSnapshot("fit2dpars_GG.root",ws->components(),kFALSE);
    
    NSigNP_static = ws->var("NSigNP")->getVal();
    NSigPR_static = ws->var("NSigPR")->getVal();
    ErrNP_static = ws->var("NSigNP")->getError();
    ErrPR_static = ws->var("NSigPR")->getError();

    Bfrac_static = NSigNP_static/(NSigNP_static + NSigPR_static);
    BfracErr_static = sqrt(pow(NSigNP_static*ErrPR_static,2) + pow(NSigPR_static*ErrNP_static,2))/pow(NSigNP_static + NSigPR_static,2);
  }

  // drawResults(ws,isGG,!prefitMass,rb2,prange,etarange);
					 
  RooRealVar tempVar1("tempVar1","tempVar1",NSigNP_static);
  RooRealVar tempVar2("tempVar2","tempVar2",NBkg_static);

  RooPlot *mframe = ws->var("JpsiMass")->frame();

  string partTit, partFile;
  if (isGG == 0) { partTit = " glb-glb "; partFile = "GG"; }
  else if (isGG == 1) { partTit = " glb-trk "; partFile = "GT"; }
  else { partTit = " all "; partFile = "ALL"; }

  titlestr = "2D fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |eta| = " + etarange;
  mframe->SetTitle(titlestr.c_str());

  reddata1->plotOn(mframe,DataError(RooAbsData::SumW2));

  if (prefitMass) {
    ws->pdf("totPDF")->plotOn(mframe,Components("expFunct"),LineColor(kBlue),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tempPDF("tempPDF","tempPDF",RooArgList(*(ws->pdf("sigCBGauss")),*(ws->pdf("expFunct"))),RooArgList(tempVar1,tempVar2));
    tempPDF.plotOn(mframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(mframe,Components("totsigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  TCanvas c1;
  c1.cd();mframe->Draw();
  titlestr = "pictures/data4Gauss/2D_" + partFile + "massfit_pT" + prange + "_eta" + etarange + ".gif";
  c1.SaveAs(titlestr.c_str());

  RooPlot *tframe = ws->var("Jpsict")->frame();

  titlestr = "2D fit for" + partTit + "muons (c  #tau projection), p_{T} = " + prange + " GeV/c and |eta| = " + etarange;
  tframe->SetTitle(titlestr.c_str());

  reddata1->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb2));

  if (prefitMass) {
    ws->pdf("totPDF")->plotOn(tframe,Components("bkgctauTOT"),LineColor(kBlue),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tempPDF2("tempPDF2","tempPDF2",RooArgList(*(ws->pdf("sigNP")),*(ws->pdf("bkgctauTOT"))),RooArgList(tempVar1,tempVar2));
    tempPDF2.plotOn(tframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(tframe,Components("totsigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(tframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  titlestr = "pictures/data4Gauss/2D_" + partFile + "timefit_pT" + prange + "_eta" + etarange + "_Lin.gif";
  c2.SaveAs(titlestr.c_str());
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  titlestr = "pictures/data4Gauss/2D_" + partFile + "timefit_pT" + prange + "_eta" + etarange + "_Log.gif";
  c2.SaveAs(titlestr.c_str());
 
  RooPlot *tframe1 = ws->var("Jpsict")->frame();

  titlestr = "2D fit for" + partTit + "muons (c  #tau projection, mass sidebands), p_{T} = " + prange + " GeV/c and |eta| = " + etarange;
  tframe->SetTitle(titlestr.c_str());

  reddataSB->plotOn(tframe1,DataError(RooAbsData::SumW2),Binning(rb2));

  ws->pdf("bkgctauTOT")->plotOn(tframe1,Normalization(reddataSB->sumEntries(),RooAbsReal::NumEvent));
  
  TCanvas c3;
  c3.cd();
  c3.cd();tframe1->Draw();
  titlestr = "pictures/data4Gauss/2D_" + partFile + "timeSide_pT" + prange + "_eta" + etarange + "_Lin.gif";
  c3.SaveAs(titlestr.c_str());
  c3.SetLogy(1);
  c3.cd();tframe1->Draw();
  titlestr = "pictures/data4Gauss/2D_" + partFile + "timeSide_pT" + prange + "_eta" + etarange + "_Log.gif";
  c3.SaveAs(titlestr.c_str()); 

  cout << endl << "J/psi yields:" << endl;
  cout << "PROMPT :     Fit : " << NSigPR_static << " +/- " << ErrPR_static << endl;
  cout << "NON-PROMPT : Fit : " << NSigNP_static << " +/- " << ErrNP_static << endl;
  cout << "B fraction : Fit : " << Bfrac_static << " +/- " << BfracErr_static << endl;

  char oFile[200];
  sprintf(oFile,"results/data4Gauss/results2D%s_pT%s_eta%s.txt",partFile.c_str(),prange.c_str(),etarange.c_str());

  ofstream outputFile(oFile);
  outputFile << "PR " << 0. << " " << NSigPR_static << " " << ErrPR_static << endl;
  outputFile << "NP " << 0. << " " << NSigNP_static << " " << ErrNP_static << endl;
  outputFile << "BF " << 0. << " " << Bfrac_static << " " << BfracErr_static << endl;
  outputFile << endl; 

  return 1;
}
