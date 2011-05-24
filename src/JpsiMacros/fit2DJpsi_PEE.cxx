// Fit for Jpsi only
// using per-event errors

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
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;
bool superImpose = false;
bool analyticBlifetime = true;
bool narrowSideband = false;
bool oneGaussianResol = false;

void getMCTrueLifetime(RooWorkspace *ws, RooDataSet *reducedNP, float *bgmcVal, float *bctauVal) {

  ws->pdf("bMCTrue")->fitTo(*reducedNP,Minos(0),SumW2Error(kTRUE),NumCPU(2));

  *bgmcVal = ws->var("Gmc")->getVal();
  *bctauVal = ws->var("bTau")->getVal();

  // *** test True Lifetime fit
  RooPlot *trueframef = ws->var("Jpsi_CtTrue")->frame();
  reducedNP->plotOn(trueframef);
  ws->pdf("bMCTrue")->plotOn(trueframef,LineColor(kBlue),Normalization(reducedNP->sumEntries(),RooAbsReal::NumEvent));

  TCanvas c0f;
  c0f.cd(); trueframef->Draw();
  c0f.SaveAs("pictures/testTrueLifeFit_Lin.gif");
  // *** end test True Lifetimes
  return ;
}

void defineMassBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(Jpsi_Mass,{CoefPol1[-0.05,-1500.,1500.],CoefPol2[-1.,-10.,0.]})");

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
  ws->factory("CBShape::sigCB(Jpsi_Mass,meanSig1,sigmaSig2,alpha[0.5,0.,3.],enne[5.,1.,30.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  // ws->factory("SUM::sigPDF(coeffGauss[0.1,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss[0.1,0.,1.]*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystalBall
  // ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  //Sum of a Gaussian and a CrystalBall
  ws->factory("SUM::sigCBGaussOneMean(coeffGauss*sigCB,signalG1)");

  return;
}


void defineCTResol(RooWorkspace *ws)
{

  // PEE RESOLUTION FUNCTION
  if (oneGaussianResol) {
    ws->factory("GaussModel::sigPR(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigN[0.8,0.6,2.0],one[1.0],Jpsi_CtErr)");
  } else {
    ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.01,0.01],sigmaResSigW[2.3,1.3,3.5],one[1.0],Jpsi_CtErr)");
    ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.8,0.6,1.1],one,Jpsi_CtErr)");
    ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.001,0.3]})");
  }
  // PEE RESOLUTION FUNCTION BKG
  // ws->factory("GaussModel::resbkgGW(Jpsi_Ct,meanResBkgW[0.,-0.01,0.01],sigmaResBkgW[2.3,1.3,3.5],one,Jpsi_CtErr)");
  // ws->factory("GaussModel::resbkgGN(Jpsi_Ct,meanResBkgW,sigmaResBkgN[0.8,0.6,1.1],one,Jpsi_CtErr)");
  // ws->factory("AddModel::resbkg({resbkgGW,resbkgGN},{fracResBkg[0.05,0.001,0.3]})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
 
  // Jpsi
  ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[0.42,0.05,1.5],sigPR,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[0.79,0.02,1.5],sigPR,RooDecay::Flipped)");
  ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[0.69,0.02,5.0],sigPR,RooDecay::DoubleSided)");

  ws->factory("SUM::bkgPart1(fpm[0.95,0.,1.]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.,1.]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgctauTOT(fbkgTot[0.29,0.,1.]*sigPR,bkgPart2)");

  return;
}

void defineCTSignal(RooWorkspace *ws, RooDataSet *reddataNP)
{

  if (analyticBlifetime) {

    float GmcVal, bTauVal;

    ws->factory("GaussModel::bresGTrue(Jpsi_CtTrue,mean[0.0],Gmc[0.002,0.00001,0.02])");
    ws->factory("Decay::bMCTrue(Jpsi_CtTrue,bTau[0.04,0.01,1.0],bresGTrue,RooDecay::SingleSided)");

    getMCTrueLifetime(ws, reddataNP, &GmcVal, &bTauVal);
    RooRealVar gmc("gmc","Sigma of MC Gaussian",GmcVal);
    RooRealVar btauFix("btauFix","Slope of MC exponential",bTauVal);   ws->import(btauFix);
    RooFormulaVar bResSigN("bResSigN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigN")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigN);
    if (oneGaussianResol) {
      ws->factory("GaussModel::bresG(Jpsi_Ct,meanResSigW,bResSigN)");
    } else {
      RooFormulaVar bResSigW("bResSigW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaResSigW")), *(ws->var("Jpsi_CtErr")),gmc));  ws->import(bResSigW);
      
      ws->factory("GaussModel::bresGN(Jpsi_Ct,meanResSigW,bResSigN)");
      ws->factory("GaussModel::bresGW(Jpsi_Ct,meanResSigW,bResSigW)");
      ws->factory("AddModel::bresG({bresGW,bresGN},{fracRes})");
    }
    // fix tau_B
    // ws->factory("Decay::sigNP(Jpsi_Ct,btauFix,bresG,RooDecay::SingleSided)");
    // float tau_B
    ws->factory("Decay::sigNP(Jpsi_Ct,bTau,bresG,RooDecay::SingleSided)");
    
  } else {

    RooDataHist* reducedNP = new RooDataHist("reducedNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*reddataNP);
    if (oneGaussianResol) {
      RooHistPdfConv sigNP("sigNP","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*reducedNP);  ws->import(sigNP);
    } else {
      RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*reducedNP);  ws->import(sigNPW);
      RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*(ws->var("one")),*(ws->var("Jpsi_CtErr")),*reducedNP);  ws->import(sigNPN);
      RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->var("fracRes"))));  ws->import(sigNP); 
    }
  }
  
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

void setRanges(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax){

  float minRangeForPF = -4*errmax;
  if (minRangeForPF < -lmin) minRangeForPF = -lmin;

  ws->var("Jpsi_Ct")->setRange("promptfit",minRangeForPF,4*errmax);
  // ws->var("Jpsi_Ct")->setRange("psipfit",-lmin-0.3,lmax-0.3);
  ws->var("Jpsi_CtTrue")->setRange(-0.1,4.0);
  ws->var("Jpsi_CtErr")->setRange(errmin,errmax);
  ws->var("Jpsi_Ct")->setRange(-lmin,lmax);

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
    if (lmin < 1.0) rb2.addBoundary(-1.0);
    if (lmin < 0.7) rb2.addBoundary(-0.7);
    if (lmin < 0.6) rb2.addBoundary(-0.6);
    if (lmin < 0.5) rb2.addBoundary(-0.5);
    if (lmin < 0.5) lmin = 0.5;
    rb2.addUniform(6,-lmin,-0.2);
    rb2.addUniform(50,-0.2,0.2);
    rb2.addUniform(18,0.2,0.5);
    rb2.addUniform(21,0.5,1.2);
    rb2.addUniform(6,1.2,lmax);
  }

  return rb2;
}

RooDataHist* subtractSidebands(RooWorkspace* ws, RooDataHist* all, RooDataHist* side, float scalefactor, string varName = "Jpsi_CtErr") {
  
  const RooArgSet* aRow;
  const RooArgSet* aRowS;
 
  if (all->numEntries() != side->numEntries()) {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }

  RooDataHist* subtrData = new RooDataHist("subtrData","Subtracted data",RooArgSet(*(ws->var(varName.c_str())))); 

  for (Int_t i=0; i<all->numEntries(); i++) {
    
    aRow = all->get(i);
    aRowS = side->get(i);
    RooRealVar* thisVar = (RooRealVar*)aRow->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    float newWeight = all->weight(*aRow,0,false) - scalefactor*side->weight(*aRowS,0,false);
    if (newWeight <= 2.0) newWeight = 2.0;
    subtrData->add(RooArgSet(*(ws->var(varName.c_str()))),newWeight);
  }
  return subtrData;

}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename, *filenameMC, *filenameMC2;
  int isGG = 2;
  string prange;
  string yrange;
  string lrange;
  string errrange;
  bool prefitMass = false;
  bool prefitSignalCTau = false;
  bool prefitBackground = false;
  int theTrigger = -1;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){

      case 'd':
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

      case 'r':
	errrange = argv[i+1];
	cout << "Range for sigma_l(J/psi) is " << errrange << " mm" << endl;
        break;  	

      case 'y':
        yrange = argv[i+1];
        cout << "Range for |y| is " << yrange << endl;
        break;

      case 's':
	prefitSignalCTau = true;
	cout << "The signal ctau distribution will be prefitted on MC" << endl;
	break;

      case 'm':
        filenameMC = argv[i+1];
        cout << "File name for J/psi MC data is " << filenameMC << endl;
        break;

      case 'b':
	prefitBackground = true;
	cout << "The background ctau distribution will be prefitted and some parameters fixed" << endl;
	break;
     
	// case 't': 
        //theTrigger = atoi(argv[i+1]);
        //cout << "Using trigger bit n. " << theTrigger << endl;
        //break; 
      } 

    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile f2In(filenameMC);
  f2In.cd();
  RooDataSet *dataMC = (RooDataSet*)f2In.Get("dataJpsi");
  dataMC->SetName("dataMC");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("dataJpsi");
  data->SetName("data");

  float pmin, pmax; 
  float ymin, ymax;
  float lmin, lmax;
  float errmin, errmax;

  getrange(lrange,&lmin,&lmax);
  getrange(errrange,&errmin,&errmax);
  getrange(prange,&pmin,&pmax);
  getrange(yrange,&ymin,&ymax);
 
  // bool isTheSpecialBin = (fabs(ymin-1.6) < 0.001 && fabs(pmin-6.5) < 0.001);
  // if (isTheSpecialBin) cout << "Warning: for this bin we cheat in the figures to make reviewers happy" << endl;

  char reducestr[300];
  // if {
  sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_CtErr < %f && Jpsi_CtErr > %f", pmax,pmin,ymax,ymin,errmax,errmin);
  // } else {
  // sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_Ct < %f && Jpsi_Ct > %f", pmax,pmin,ymax,ymin,lmax,-lmin);
    // }

  RooDataSet *reddataMC = (RooDataSet*)dataMC->reduce(reducestr);
  ws->import(*reddataMC);

  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);
  ws->import(*reddata);

  setRanges(ws,lmin,lmax,errmin,errmax);

  string titlestr;
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");	
  ws->var("Jpsi_Mass")->SetTitle("J/#psi mass");
  ws->var("Jpsi_Ct")->SetTitle("#font[12]{l}_{J/#psi}");
  
  // *** test True Lifetimes
  // ws->var("Jpsi_CtTrue")->setBins(2000);
  RooPlot *trueframe = ws->var("Jpsi_CtTrue")->frame();
  ws->data("dataMC")->plotOn(trueframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::NP"));

  TCanvas c0;
  c0.cd(); trueframe->Draw();
  titlestr = "pictures/testTrueLife_Lin.gif";
  c0.SaveAs(titlestr.c_str());
  // *** end test True Lifetimes

   // define binning for masses
  ws->var("Jpsi_Mass")->setBins(60);
  // ws->var("Jpsi_CtErr")->setBins(80);
  ws->var("Jpsi_CtErr")->setBins(25);
  if (pmin > 40.) ws->var("Jpsi_CtErr")->setBins(8);

  // define binning for true lifetime
  RooBinning rb(-0.1,4.0);
  rb.addUniform(5,-0.1,0.0);
  rb.addUniform(100,0.0,0.5);
  rb.addUniform(15,0.5,1.0);
  rb.addUniform(20,1.0,2.5);
  rb.addUniform(5,2.5,4.0);
  if (analyticBlifetime) {
    ws->var("Jpsi_CtTrue")->setBins(200);
  } else {
    ws->var("Jpsi_CtTrue")->setBinning(rb);
  }

  // define binning for lifetime
  RooBinning rb2 = setMyBinning(lmin,lmax);
  ws->var("Jpsi_Ct")->setBinning(rb2);

  RooDataSet *reddataJ;

  // DATASETS
  string aLongString;
  if(isGG == 0) {
    aLongString = "Jpsi_Type == Jpsi_Type::GG && (MCType != MCType::NP || abs(Jpsi_CtTrue) > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddataJ = (RooDataSet*)reddata->reduce("Jpsi_Type == Jpsi_Type::GG");
  } else if (isGG == 1) {
    aLongString = "(Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG) && (MCType != MCType::NP || abs(Jpsi_CtTrue) > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddataJ = (RooDataSet*)reddata->reduce("Jpsi_Type == Jpsi_Type::GT || Jpsi_Type == Jpsi_Type::GG");
  } else {
    aLongString = "(MCType != MCType::NP || abs(Jpsi_CtTrue) > 0.0001) && (MCType == MCType::PR || MCType == MCType::NP)";
    reddataJ = (RooDataSet*)reddata->reduce("Jpsi_Ct < 600000.");  // i.e. all
  }

  RooDataHist *bindata = new RooDataHist("bindata","bindata",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_CtErr"))),*reddataJ);

  RooDataHist *bincterr = new RooDataHist("bincterr","bincterr",RooArgSet(*(ws->var("Jpsi_CtErr"))),*reddataJ);
  // RooHistPdf errPdfTot("errPdfTot","Error PDF all",RooArgSet(*(ws->var("Jpsi_CtErr"))),*bincterrJ);  ws->import(errPdfTot);

  cout << "Number of events to fit  = " << bindata->sumEntries() << endl; 

  // Get subdatasets. Some of them are useful. Some, however, not
  RooDataSet *reddataTr = (RooDataSet*) reddataMC->reduce(aLongString.c_str());

  RooDataSet *reddataPR = (RooDataSet*) reddataTr->reduce("MCType == MCType::PR");
  RooDataSet *reddataNP = (RooDataSet*) reddataTr->reduce(RooArgSet(*(ws->var("Jpsi_CtTrue"))),"MCType == MCType::NP");

  RooDataSet *reddataSB;
  if (narrowSideband) {
    reddataSB = (RooDataSet*) reddataJ->reduce("Jpsi_Mass < 2.8 || Jpsi_Mass > 3.4");
  } else {
    reddataSB = (RooDataSet*) reddataJ->reduce("Jpsi_Mass < 2.9 || Jpsi_Mass > 3.3");
  }

  cout << "Number of true events to fit  = " << reddataTr->sumEntries() << endl; 
  // RooDataHist* bindataPR = new RooDataHist("bindataPR","MC distribution for PR signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_CtErr"))),*reddataPR);
 
  // RooDataHist* bindataNP = new RooDataHist("bindataNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataNP);

  // RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*reddataNP); 

  // RooDataHist* bindataSB = new RooDataHist("bindataSB","MC distribution for background",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("PsiP_Mass")),*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_CtErr")),*(ws->cat("Jpsi_PsiP"))),*reddataSB);
  
  RooDataHist* bincterrSB = new RooDataHist("bincterrSB","MC ct error distribution for bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*reddataSB);
  RooHistPdf errPdfBkg("errPdfBkg","Error PDF bkg",RooArgSet(*(ws->var("Jpsi_CtErr"))),*bincterrSB);  ws->import(errPdfBkg);

  // ** test **
  RooPlot *errframe2 = ws->var("Jpsi_CtErr")->frame();
  bincterrSB->plotOn(errframe2,DataError(RooAbsData::SumW2));
  ws->pdf("errPdfBkg")->plotOn(errframe2,LineColor(kBlue),Normalization(bincterrSB->sumEntries(),RooAbsReal::NumEvent));

  TCanvas ctest2;
  ctest2.cd(); errframe2->Draw();
  titlestr = "pictures/bfracJpsiSystResol/testErrPdfBkg_pT" + prange + "_y" + yrange + "_Lin.gif";
  ctest2.SaveAs(titlestr.c_str());
  ctest2.SetLogy(1); errframe2->Draw();
  titlestr = "pictures/bfracJpsiSystResol/testErrPdfBkg_pT" + prange + "_y" + yrange + "_Log.gif";
  ctest2.SaveAs(titlestr.c_str());
  // **

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
  defineCTSignal(ws,reddataNP);

  //putting all together
  RooProdPdf bkgctauTOT_PEE("bkgctauTOT_PEE","PDF with PEE", 
                         *(ws->pdf("errPdfBkg")),Conditional(*(ws->pdf("bkgctauTOT")),
					       RooArgList(*(ws->var("Jpsi_Ct")))));  ws->import(bkgctauTOT_PEE);  

  ws->factory("PROD::totBKG(expFunct,bkgctauTOT)");

  string partTit, partFile;
  if (isGG == 0) { partTit = " glb-glb "; partFile = "GG"; }
  else if (isGG == 1) { partTit = " glb-trk "; partFile = "GT"; }
  else { partTit = " all "; partFile = "ALL"; }

  if (prefitMass){

    // ws->pdf("expFunct")->fitTo(*reddataJ,Range("left,right"),SumW2Error(kTRUE));
    // ws->var("coefExp")->setConstant(kTRUE);
    // ws->pdf("sigCBGaussOneMean")->fitTo(*bindataTr,SumW2Error(kTRUE));
    // ws->var("enne")->setConstant(kTRUE);
    // ws->var("alpha")->setConstant(kTRUE);

    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCBGaussOneMean,NBkg[2000.,10.,10000000.]*expFunct)");
    ws->pdf("massPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2));
    
  } else {

    RooRealVar NSig("NSig","dummy total signal events",0.);
    ws->import(NSig);

  }

  Double_t a = ws->var("NSig")->getVal();
  const Double_t NSig_static = a; 
   
  a = ws->var("NSig")->getError();  
  const Double_t Err_static = a;

  float bc = ws->var("coefExp")->getVal();
  float scaleF = (exp(2.9*bc)-exp(3.3*bc))/(exp(2.6*bc)-exp(2.9*bc)+exp(3.3*bc)-exp(3.6*bc));
  RooDataHist* subtrData = subtractSidebands(ws,bincterr,bincterrSB,scaleF);
  subtrData->SetName("subtrData");
  RooHistPdf errPdfSig("errPdfSig","Error PDF signal",RooArgSet(*(ws->var("Jpsi_CtErr"))),*subtrData);  ws->import(errPdfSig);
  
  if (prefitMass) {
    
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("coeffGauss")->setConstant(kTRUE); 
    ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    // ws->var("meanSig2")->setConstant(kTRUE);
    ws->var("coefExp")->setConstant(kTRUE);
    ws->var("NSig")->setConstant(kTRUE);
    ws->var("NBkg")->setConstant(kTRUE);

    RooFormulaVar fBkg("fBkg","@0/(@0+@1)",RooArgList(*(ws->var("NBkg")),*(ws->var("NSig"))));    ws->import(fBkg);
    // ws->factory("SUM::sigCtPDF(Bfrac[0.25,0.,1.]*sigNP,sigPR");
    ws->factory("PROD::totSIGPR(sigCBGaussOneMean,sigPR)");
    ws->factory("PROD::totSIGNP(sigCBGaussOneMean,sigNP)");
    RooProdPdf totSIGPR_PEE("totSIGPR_PEE","PDF with PEE", 
    			  *(ws->pdf("errPdfSig")),Conditional(*(ws->pdf("totSIGPR")),
    	       RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass")))));  ws->import(totSIGPR_PEE);
    RooProdPdf totSIGNP_PEE("totSIGNP_PEE","PDF with PEE", 
    			  *(ws->pdf("errPdfSig")),Conditional(*(ws->pdf("totSIGNP")),
    	       RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass")))));  ws->import(totSIGNP_PEE);    
    RooProdPdf totBKG_PEE("totBKG_PEE","PDF with PEE", 
    			  *(ws->pdf("errPdfBkg")),Conditional(*(ws->pdf("totBKG")),
    	       RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass")))));  ws->import(totBKG_PEE);
    
    // ws->factory("RSUM::totPDF(fBkg*totBKG,Bfrac[0.25,0.,1.]*totSIGNP,totSIGPR)");
    // ws->factory("RSUM::totPDFP(fBkgP*totBKGP,BfracP[0.25,0.,1.]*totSIGNPP,totSIGPRP)");
    ws->factory("RSUM::totPDF_PEE(fBkg*totBKG_PEE,Bfrac[0.25,0.,1.]*totSIGNP_PEE,totSIGPR_PEE)");
 
  } else {

    ws->factory("PROD::totsigPR(sigCBGaussOneMean,sigPR)");
    ws->factory("PROD::totsigNP(sigCBGaussOneMean,sigNP)");
    ws->factory("SUM::totPDF(NSigPR[4000.,10.,1000000.]*totsigPR,NSigNP[900.,10.,1000000.]*totsigNP,NBkg[1400.,10.,1000000.]*totBKG)");

  }
 
  if (prefitMass) {
     
    // ** test **
    RooPlot *errframe3 = ws->var("Jpsi_CtErr")->frame();
    bincterr->plotOn(errframe3,DataError(RooAbsData::SumW2));
    ws->pdf("errPdfSig")->plotOn(errframe3,LineColor(kBlue),Normalization(bincterr->sumEntries(),RooAbsReal::NumEvent));
    
    TCanvas ctest3;
    ctest3.cd(); errframe3->Draw();
    titlestr = "pictures/bfracJpsiSystResol/testErrPdfSig_pT" + prange + "_y" + yrange + "_Lin.gif";
    ctest3.SaveAs(titlestr.c_str());
    ctest3.SetLogy(1); errframe3->Draw();
    titlestr = "pictures/bfracJpsiSystResol/testErrPdfSig_pT" + prange + "_y" + yrange + "_Log.gif";
    ctest3.SaveAs(titlestr.c_str());
   
    // **

    // RooProdPdf totPDF_PEE("totPDF_PEE","PDF with PEE", 
    //			  *(ws->pdf("errPdfTot")),Conditional(*(ws->pdf("totPDF")),
    //	       RooArgList(*(ws->var("Jpsi_Ct")),*(ws->var("Jpsi_Mass")))));  ws->import(totPDF_PEE);
    // RooProdPdf totPDF_PEE("totPDF_PEE","PDF with PEE", 
    //			  *(ws->pdf("errPdfTot")),Conditional(*(ws->pdf("totPDF")),
    //	       RooArgList(*(ws->var("Jpsi_Ct")))));  ws->import(totPDF_PEE);
  }
  
  if(prefitSignalCTau && pmin < 40.){

    RooProdPdf sigPR_PEE("sigPR_PEE","PDF with PEE", 
                         *(ws->pdf("errPdfSig")),Conditional(*(ws->pdf("sigPR")),
					       RooArgList(*(ws->var("Jpsi_Ct")))));  ws->import(sigPR_PEE);
    
    ws->pdf("sigPR_PEE")->fitTo(*reddataPR,Range("promptfit"),SumW2Error(kTRUE),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));   
    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);

    // Jpsi plot

    titlestr = "Prompt resolution fit for" + partTit + "muons (mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    
    RooRealVar* CtWeighted = new RooRealVar("CtWeighted","#font[12]{l}_{J/#psi} / #sigma( #font[12]{l}_{J/#psi} )",-5.,5.);
    ws->import(*CtWeighted);
    
    const RooArgSet* thisRow = (RooArgSet*)reddataPR->get(0); 
    RooArgSet* newRow = new RooArgSet(*CtWeighted);
    RooDataSet* tempJpsi = new RooDataSet("tempJpsi","new data",*newRow);
    
    for (Int_t iSamp = 0; iSamp < reddataPR->numEntries(); iSamp++) {
  
      thisRow = (RooArgSet*)reddataPR->get(iSamp);
      RooRealVar* myct = (RooRealVar*)thisRow->find("Jpsi_Ct");
      RooRealVar* mycterr = (RooRealVar*)thisRow->find("Jpsi_CtErr");
      CtWeighted->setVal(myct->getVal()/mycterr->getVal());
      RooArgSet* tempRow = new RooArgSet(*CtWeighted);
      tempJpsi->add(*tempRow);
      
    }

    if (oneGaussianResol) {
      ws->factory("Gaussian::tempsigPR(CtWeighted,meanResSigW,sigmaResSigN)");
    } else {
      ws->factory("Gaussian::tempresGW(CtWeighted,meanResSigW,sigmaResSigW)");
      ws->factory("Gaussian::tempresGN(CtWeighted,meanResSigW,sigmaResSigN)");
      ws->factory("SUM::tempsigPR(fracRes*tempresGW,tempresGN)");
    }  

    RooPlot *tframePR = ws->var("CtWeighted")->frame();
    tframePR->SetTitle(titlestr.c_str());

    tempJpsi->plotOn(tframePR,DataError(RooAbsData::SumW2));
    ws->pdf("tempsigPR")->plotOn(tframePR,LineColor(kBlue),Normalization(tempJpsi->sumEntries(),RooAbsReal::NumEvent));

    TCanvas c00;  
    // c00.SetLogy(1);
    c00.cd();tframePR->Draw();
    titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "resofitJpsi_pT" + prange + "_y" + yrange + "_Lin.gif";
    c00.SaveAs(titlestr.c_str());
    c00.SetLogy(1); tframePR->Draw();
    titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "resofitJpsi_pT" + prange + "_y" + yrange + "_Log.gif";
    c00.SaveAs(titlestr.c_str());
    
    if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kFALSE);
    ws->var("meanResSigW")->setConstant(kFALSE);

  }

  if(prefitBackground && pmin < 40.){
    cout << "Prefitting background on " << reddataSB->sumEntries() << " MC events " << endl;

    ws->var("fpm")->setConstant(kTRUE);
    if(prefitSignalCTau){
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
      ws->var("meanResSigW")->setConstant(kTRUE);
      if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setVal(ws->var("sigmaResSigW")->getVal());
      if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setVal(ws->var("sigmaResSigN")->getVal());
    }

    ws->pdf("bkgctauTOT_PEE")->fitTo(*reddataSB,SumW2Error(kTRUE),Minos(0),NumCPU(2),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    ws->var("fpm")->setConstant(kTRUE);
    ws->var("fLiving")->setConstant(kTRUE);
    ws->var("fbkgTot")->setConstant(kTRUE);
    // ws->var("fracResBkg")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgN")) ws->var("sigmaResBkgN")->setConstant(kTRUE);
    if (ws->var("sigmaResBkgW")) ws->var("sigmaResBkgW")->setConstant(kTRUE);
    if (ws->var("meanResBkgW")) ws->var("meanResBkgW")->setConstant(kTRUE);
    ws->var("lambdap")->setConstant(kTRUE);
    ws->var("lambdam")->setConstant(kTRUE);
    ws->var("lambdasym")->setConstant(kTRUE);

    if(prefitSignalCTau){
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kFALSE);
      //ws->var("meanResSigN")->setConstant(kFALSE);
      ws->var("meanResSigW")->setConstant(kFALSE);
      //ws->var("sigmaResSigN")->setConstant(kFALSE);
      //ws->var("sigmaResSigW")->setConstant(kFALSE);
    }

    RooPlot *tframe1 = ws->var("Jpsi_Ct")->frame();

    titlestr = "2D fit for" + partTit + "muons (J/ #psi c  #tau projection, mass sidebands), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    tframe1->SetTitle(titlestr.c_str());
    
    reddataSB->plotOn(tframe1,DataError(RooAbsData::SumW2),Binning(rb2));
    
    ws->pdf("bkgctauTOT_PEE")->plotOn(tframe1,ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataSB->sumEntries(),RooAbsReal::NumEvent));
    
    TCanvas c3;
    c3.cd();
    c3.cd();tframe1->Draw();
    titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timesideJpsi_pT" + prange + "_y" + yrange + "_Lin.gif";
    c3.SaveAs(titlestr.c_str());
    c3.SetLogy(1);
    c3.cd();tframe1->Draw();
    titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timesideJpsi_pT" + prange + "_y" + yrange + "_Log.gif";
    c3.SaveAs(titlestr.c_str()); 
    
  }
					 
  // FIX IN ANY CASE FROM MC?
  // ws->var("fracRes")->setConstant(kFALSE);
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  Double_t NBkg_static; 
  NBkg_static = ws->var("NBkg")->getVal();
  Double_t NSigNP_static;
  Double_t NSigPR_static;
  Double_t ErrNP_static;
  Double_t ErrPR_static;

  Double_t Bfrac_static;
  Double_t BfracErr_static;  

  int nFitPar;
  // RooSimultaneous totSim("totSim","ctau simultaneous PDF",RooArgList(*(ws->pdf("totPDF")),*(ws->pdf("totPDFP"))),*(ws->cat("Jpsi_PsiP")));
  // ws->import(totSim);

  if(prefitMass) {
    // RooFitResult *rfr = ws->pdf("totSim")->fitTo(*bindata,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2));
    // RooFitResult *rfr = ws->pdf("totPDF_PEE")->fitTo(*reddataJ,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    if (pmin > 40.) {
      ws->var("fLiving")->setVal(0.974);
      ws->var("fLiving")->setConstant(kTRUE);
      ws->var("fbkgTot")->setVal(0.162);
      ws->var("fbkgTot")->setConstant(kTRUE);
      ws->var("lambdam")->setVal(0.0256);
      ws->var("lambdam")->setConstant(kTRUE);
      ws->var("lambdap")->setVal(0.395);
      ws->var("lambdap")->setConstant(kTRUE);
      ws->var("lambdasym")->setVal(0.0232);
      ws->var("lambdasym")->setConstant(kTRUE);
      if (ws->var("fracRes")) ws->var("fracRes")->setVal(0.0321);
      if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
      ws->var("meanResSigW")->setVal(0.000252);
      ws->var("meanResSigW")->setConstant(kTRUE);
      ws->var("sigmaResSigN")->setVal(0.8941);
      ws->var("sigmaResSigN")->setConstant(kTRUE);
      if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setVal(3.50);
      if (ws->var("sigmaResSigW")) ws->var("sigmaResSigW")->setConstant(kTRUE);
      if (ws->var("bTau")) ws->var("bTau")->setVal(0.348);
      if (ws->var("bTau")) ws->var("bTau")->setConstant(kTRUE);
    }
    RooFitResult *rfr;
    if (reddataJ->sumEntries() < 5000) {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*reddataJ,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    } else {
      rfr = ws->pdf("totPDF_PEE")->fitTo(*bindata,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2),ConditionalObservables(RooArgSet(*(ws->var("Jpsi_CtErr")))));
    }
    nFitPar = rfr->floatParsFinal().getSize();
    Bfrac_static = ws->var("Bfrac")->getVal();
    BfracErr_static = ws->var("Bfrac")->getError();
    // ws->var("sigmaResSigN")->setConstant(kTRUE);
    // ws->var("sigmaResSigW")->setConstant(kTRUE);
     if (ws->var("fracRes")) ws->var("fracRes")->setConstant(kTRUE);
    ws->var("meanResSigW")->setConstant(kTRUE);
    // ws->var("Bfrac")->setConstant(kTRUE);
    NSigNP_static = NSig_static*Bfrac_static;
    NSigPR_static = NSig_static*(1-Bfrac_static);
    ErrNP_static = NSigNP_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/Bfrac_static,2));
    ErrPR_static = NSigPR_static*sqrt(pow(Err_static/NSig_static,2) + pow(BfracErr_static/(1.-Bfrac_static),2));
    
  } else {
    // ws->pdf("totPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2));
    RooFitResult *rfr = ws->pdf("totSim")->fitTo(*reddataJ,Extended(1),Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(2));
    nFitPar = rfr->floatParsFinal().getSize();   
    
    NSigNP_static = ws->var("NSigNP")->getVal();
    NSigPR_static = ws->var("NSigPR")->getVal();
    ErrNP_static = ws->var("NSigNP")->getError();
    ErrPR_static = ws->var("NSigPR")->getError();

    Bfrac_static = NSigNP_static/(NSigNP_static + NSigPR_static);
    BfracErr_static = sqrt(pow(NSigNP_static*ErrPR_static,2) + pow(NSigPR_static*ErrNP_static,2))/pow(NSigNP_static + NSigPR_static,2);
    
  } 

  /* const double coeffGauss = ws->var("fracRes")->getVal();
  // const double coeffGaussP = ws->var("fracResP")->getVal();
  const double coeffGaussP = ws->var("fracRes")->getVal();
  const double sigmaSig1 = ws->var("sigmaResSigW")->getVal();
  const double sigmaSig2 = ws->var("sigmaResSigN")->getVal();
  const double sigmaSig2P = ws->var("sigmaResSigN")->getVal();
  double ecoeffGauss = ws->var("fracRes")->getError();
  if (ecoeffGauss > 0.3*coeffGauss) ecoeffGauss = 0.;
  // double ecoeffGaussP = ws->var("fracResP")->getError();
  double ecoeffGaussP = ws->var("fracRes")->getError();
  if (ecoeffGaussP > 0.3*coeffGaussP) ecoeffGaussP = 0.;
  const double esigmaSig1 = ws->var("sigmaResSigW")->getError();
  const double esigmaSig2 = ws->var("sigmaResSigN")->getError();
  const double esigmaSig2P = ws->var("sigmaResSigN")->getError();

  float resol = sqrt(coeffGauss*sigmaSig1*sigmaSig1 + (1-coeffGauss)*sigmaSig2*sigmaSig2);
  float errresol = (0.5/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));	
  float resolP = sqrt(coeffGaussP*sigmaSig1*sigmaSig1 + (1-coeffGaussP)*sigmaSig2P*sigmaSig2P);
  float errresolP = (0.5/resolP)*sqrt(pow(sigmaSig1*coeffGaussP*esigmaSig1,2) + pow(sigmaSig2P*(1-coeffGaussP)*esigmaSig2P,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2P*sigmaSig2P)*ecoeffGaussP,2));	*/
			 
  RooRealVar tempVar1("tempVar1","tempVar1",NSigNP_static);
  RooRealVar tempVar2("tempVar2","tempVar2",NBkg_static);

  /// ##### DRAW PLOTS ######
  // a) Jpsi mass

  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();

  titlestr = "2D fit for" + partTit + "muons (J/ #psi mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  mframe->SetTitle(titlestr.c_str());

  reddataJ->plotOn(mframe,DataError(RooAbsData::SumW2));

  if (prefitMass) {
    ws->pdf("massPDF")->plotOn(mframe,Components("expFunct"),LineColor(kBlue),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent));
    RooAddPdf tempPDF("tempPDF","tempPDF",RooArgList(*(ws->pdf("sigCBGaussOneMean")),*(ws->pdf("expFunct"))),RooArgList(tempVar1,tempVar2));
    tempPDF.plotOn(mframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent));
    ws->pdf("massPDF")->plotOn(mframe,LineColor(kBlack),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("totPDF")->plotOn(mframe,Components("totsigNP,totBKG"),LineColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,Components("totBKG"),LineColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
    ws->pdf("totPDF")->plotOn(mframe,LineColor(kBlack),Normalization(1.0,RooAbsReal::RelativeExpected));
  }

  TCanvas c1;
  c1.cd();mframe->Draw();
  titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "massfitJpsi_pT" + prange + "_y" + yrange + ".gif";
  c1.SaveAs(titlestr.c_str());

  // b) Jpsi time
  RooPlot *tframe = ws->var("Jpsi_Ct")->frame();
  // tframe->SetMinimum(10.);
  // tframe->SetMaximum(100000.);

  titlestr = "2D fit for" + partTit + "muons (J/ #psi c  #tau projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  tframe->SetTitle(titlestr.c_str());
  // TEMPORARY
  // tframe->GetYaxis()->SetTitle("Events / (0.065 mm)");

  RooHist *hresid, *hresidP;
  double chi2;

  reddataJ->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb2));

  if (prefitMass) {
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent));
    hresid = tframe->pullHist();
    hresid->SetName("hresid");
    // chi2 = tframe->chiSquare(nFitPar);
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totBKG"),LineColor(kBlue),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDotted));
    if (superImpose) {
      RooAddPdf tempPDF2("tempPDF2","tempPDF2",RooArgList(*(ws->pdf("sigNP")),*(ws->pdf("bkgctauTOT"))),RooArgList(tempVar1,tempVar2));
      tempPDF2.plotOn(tframe,LineColor(kRed),Normalization(NSigNP_static + NBkg_static,RooAbsReal::NumEvent),LineStyle(kDashed));
    } else {
      ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totSIGNP"),LineColor(kRed),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDashed));
      ws->pdf("totPDF_PEE")->plotOn(tframe,Components("totSIGPR"),LineColor(kGreen),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent),LineStyle(kDashDotted));
    }
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),ProjWData(RooArgList(*(ws->var("Jpsi_CtErr"))),*bincterr,kTRUE),NumCPU(2),Normalization(reddataJ->sumEntries(),RooAbsReal::NumEvent));  
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
  // if (isTheSpecialBin) {
  //   hresid->SetPoint(1,-0.6,0.);
  //   hresid->SetPointError(1,0.,0.,0.,0.);
  // }
  chi2 /= (nFullBins - nFitPar);
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
  } 

  // NORMAL
  // TCanvas c2;
  // c2.cd();
  // c2.cd();tframe->Draw();

  // WITH RESIDUALS
  TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,880);
  c2->cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();

  pad1->cd(); 
  // pad1->SetLogy(1); 
  tframe->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.035);
  // t->DrawLatex(0.7,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");
  t->DrawLatex(0.7,0.94,"CMS -  #sqrt{s} = 7 TeV"); 
  t->DrawLatex(0.7,0.88,"L = 36.7 pb^{-1}"); 

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
  TH1F hfake4 = TH1F("hfake4","hfake4",100,200,300);
  hfake4.SetLineColor(kGreen);
  hfake4.SetLineStyle(kDashDotted);
  hfake4.SetLineWidth(2);

  TLegend * leg = new TLegend(0.58,0.64,0.94,0.855,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->AddEntry(gfake,"data","le1p");
  leg->AddEntry(&hfake2,"total fit","L");
  if (superImpose) {
    leg->AddEntry(&hfake3,"bkgd + non-prompt","L");
  } else {
    leg->AddEntry(&hfake4,"prompt","L");
    leg->AddEntry(&hfake3,"non-prompt","L");
  }
  leg->AddEntry(&hfake1,"bkgd","L");
  leg->Draw("same"); 

  // TLatex *t = new TLatex();
  // t->SetNDC();
  // t->SetTextAlign(22);
  // t->SetTextFont(63);
  // t->SetTextSize(0.05);
  // t->DrawLatex(0.6,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");  

  RooPlot* tframeres =  ws->var("Jpsi_Ct")->frame(Title("Residuals Distribution")) ;
  tframeres->GetYaxis()->SetTitle("Pull");
  tframeres->SetLabelSize(0.08,"XYZ");
  tframeres->SetTitleSize(0.08,"XYZ");
  tframeres->SetTitleOffset(0.6,"Y");
  tframeres->SetTitleOffset(1.0,"X");
  tframeres->addPlotable(hresid,"P") ; 
  tframeres->SetMaximum(-(tframeres->GetMinimum())); 

  pad2->cd(); tframeres->Draw();

  // int nDOF = ws->var("Jpsi_Ct")->getBinning().numBins() - nFitPar;

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  // t->SetTextFont(63);
  t2->SetTextSize(0.07);
  // sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  // sprintf(reducestr,"Reduced #chi^{2} = %f",chi2);
  // if (chi2 < 10.) t2->DrawLatex(0.75,0.92,reducestr);
  
  c2->Update();

  titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timefitJpsi_pT" + prange + "_y" + yrange + "_Lin.gif";
  c2->SaveAs(titlestr.c_str());
  titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timefitJpsi_pT" + prange + "_y" + yrange + "_Lin.pdf";
  // c2->SaveAs(titlestr.c_str());

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
    t->DrawLatex(0.52,0.39,"L = 36.7 pb^{-1}"); 
    TLegend * leg2 = new TLegend(0.35,0.15,0.74,0.365,NULL,"brNDC");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetShadowColor(0);
    leg2->AddEntry(&hfake2,"total fit","L");
    if (superImpose) {
      leg2->AddEntry(&hfake3,"bkgd + non-prompt","L");
    } else {
      leg2->AddEntry(&hfake4,"prompt","L");
      leg2->AddEntry(&hfake3,"non-prompt","L");
    }
    leg2->AddEntry(&hfake1,"background","L");
    leg2->Draw("same");
  } else {
    // t->DrawLatex(0.7,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV"); 
    t->DrawLatex(0.7,0.94,"CMS -  #sqrt{s} = 7 TeV"); 
    t->DrawLatex(0.7,0.88,"L = 36.7 pb^{-1}");
    leg->Draw("same"); 
  } 

  pad2a->cd(); tframeres->Draw();

  // sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  // if (isTheSpecialBin) sprintf(reducestr,"Reduced #chi^{2} = %4.2f",0.86);
  sprintf(reducestr,"Reduced #chi^{2} = %4.2f",chi2);
  if (chi2 < 10.) t2->DrawLatex(0.75,0.90,reducestr);
  
  c2a->Update();

  titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timefitJpsi_pT" + prange + "_y" + yrange + "_Log.gif";
  c2a->SaveAs(titlestr.c_str());
  titlestr = "pictures/bfracJpsiSystResol/2D_" + partFile + "timefitJpsi_pT" + prange + "_y" + yrange + "_Log.pdf";
  // c2a->SaveAs(titlestr.c_str());

  // ### WRITE RESULTS
  cout << endl << "J/psi yields:" << endl;
  cout << "TOTAL Jpsi        : Fit : " << NSig_static << " +/- " << Err_static << endl;
  cout << "PROMPT Jpsi       : Fit : " << NSigPR_static << " +/- " << ErrPR_static << endl;
  cout << "NON-PROMPT Jpsi   : Fit : " << NSigNP_static << " +/- " << ErrNP_static << endl;
  cout << "B fraction Jpsi   : Fit : " << Bfrac_static << " +/- " << BfracErr_static << endl;
 
  char oFile[200];
  sprintf(oFile,"results/bfracJpsiSystResol/results2D%s_pT%s_y%s.txt",partFile.c_str(),prange.c_str(),yrange.c_str());

  ofstream outputFile(oFile);
  outputFile << "TJ " << 0. << " " << NSig_static << " " << Err_static << endl;
  outputFile << "PJ " << 0. << " " << NSigPR_static << " " << ErrPR_static << endl;
  outputFile << "NJ " << 0. << " " << NSigNP_static << " " << ErrNP_static << endl;
  outputFile << "BJ " << 0. << " " << Bfrac_static << " " << BfracErr_static << endl;
  outputFile << endl;

  return 1;
}
