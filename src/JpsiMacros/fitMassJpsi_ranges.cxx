// fit for Jpsi only

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
#include "RooSimultaneous.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;
bool superImpose = false;
bool verbose = true;
int nVtx = 0;
bool cutOnVert = false;

void defineMassBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(Jpsi_Mass,{CoefPol1[-0.05,-10.,10.]})");
  ws->factory("Polynomial::CPolFunctP(Jpsi_Mass,{CoefPol1P[-0.05,-1500.,1500.],CoefPol2P[-1.,-10.,0.]})");
  // ws->factory("SUM::totMassBkgPol(sumPol[0.1,0.,1.]*CPolFunctP,CPolFunct)");

  //Exponential
  ws->factory("Exponential::expFunct(Jpsi_Mass,coefExp[-1.,-3.,1.])");
  ws->factory("Exponential::expFunctP(Jpsi_Mass,coefExpP[-1.,-3.,3.])");
  ws->factory("SUM::totMassBkgExp(sumExp[0.1,0.,1.]*expFunctP,expFunct)");

  return;
}

void defineMassSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.02,0.008,0.2])");
  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");

  // Fix Jpsi-psi' mass difference
  RooFormulaVar meanSig1P("meanSig1P","@0+0.58917",RooArgList(*(ws->var("meanSig1"))));  ws->import(meanSig1P);
  // Fix resolution scale: sigma_MJpsi/MJpsi = sigma_Mpsi'/Mpsi'
  RooFormulaVar sigmaSig1P("sigmaSig1P","@0*1.1902",RooArgList(*(ws->var("sigmaSig1"))));  ws->import(sigmaSig1P);
  RooFormulaVar sigmaSig2P("sigmaSig2P","@0*1.1902",RooArgList(*(ws->var("sigmaSig2"))));  ws->import(sigmaSig2P);
  ws->factory("Gaussian::signalG1P(Jpsi_Mass,meanSig1P,sigmaSig1P)");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(Jpsi_Mass,meanSig1,sigmaSig2)");
  ws->factory("Gaussian::signalG2OneMeanP(Jpsi_Mass,meanSig1P,sigmaSig2P)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(Jpsi_Mass,meanSig1,sigmaSig2,alpha[1.7,0.5,3.],enne[2.,1.5,10.])");
  // Same tail parameters!
  ws->factory("CBShape::sigCBP(Jpsi_Mass,meanSig1P,sigmaSig2P,alpha,enne");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  // ws->factory("SUM::sigPDF(coeffGauss[0.1,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss[0.1,0.,1.]*signalG1,signalG2OneMean)");
  ws->factory("SUM::sigPDFOneMeanP(coeffGaussP[0.1,0.,1.]*signalG1P,signalG2OneMeanP)");

  //Sum of a Gaussian and a CrystalBall
  // ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  //Sum of a Gaussian and a CrystalBall
  ws->factory("SUM::sigCBGaussOneMean(coeffGauss*sigCB,signalG1)");
  ws->factory("SUM::sigCBGaussOneMeanP(coeffGaussP*sigCBP,signalG1P)");

  return;
}


/* void defineCTResol(RooWorkspace *ws)
{

  // ONE RESOLUTION FUNCTION
  // ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.07,0.03,0.6])");
  ws->factory("GaussModel::resGW(Jpsi_Ct,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.07,0.01,0.9])");
  // ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.04,0.01,0.18])");
  ws->factory("GaussModel::resGNP(Jpsi_Ct,meanResSigWP[0.,-0.1,0.1],sigmaResSigNP[0.04,0.01,0.3])");
  ws->factory("GaussModel::resGN(Jpsi_Ct,meanResSigW,sigmaResSigN[0.04,0.01,0.3])");
  ws->factory("GaussModel::resGO(Jpsi_Ct,meanResSigW,sigmaResSigO[0.2,0.05,0.7])");
  ws->factory("GaussModel::resGM(Jpsi_Ct,meanResSigW,sigmaResSigM[0.4,0.06,0.9])");
  // ws->factory("AddModel::sigPR({resGW,resGN},{fracRes[0.05,0.,0.5]})");
  // ws->factory("AddModel::sigPR({resGW,resGO,resGM,resGN},{fracRes[0.2,0.01,0.5],fracRes2[0.02,0.0,0.30],fracRes3[0.1,0.001,0.5]})");
  ws->factory("AddModel::sigPR({resGW,resGO,resGM,resGN},{fracRes[0.2,0.01,0.9],fracRes2[0.02,0.0,0.25],fracRes3[0.1,0.0,0.15]})");
  // ws->factory("AddModel::sigPRP({resGW,resGO,resGM,resGNP},{fracResP[0.2,0.01,0.9],fracRes2,fracRes3})");
  ws->factory("AddModel::sigPRP({resGW,resGO,resGM,resGNP},{fracRes,fracRes2,fracRes3})");

  // ANOTHER RESOLUTION FUNCTION
  ws->factory("GaussModel::resbkgGW(Jpsi_Ct,meanResBkgW[0.,-0.1,0.1],sigmaResBkgW[0.05,0.005,0.5])");
  ws->factory("GaussModel::resbkgGN(Jpsi_Ct,meanResBkgW,sigmaResBkgN[0.01,0.005,0.5])");
  ws->factory("AddModel::resbkg({resbkgGW,resGO,resGM,resbkgGN},{fracResBkg[0.05,0.,1.],fracRes2,fracRes3})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
 
  // Jpsi
  ws->factory("Decay::bkg2(Jpsi_Ct,lambdap[0.42,0.05,1.5],resbkg,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3(Jpsi_Ct,lambdam[0.79,0.02,1.5],resbkg,RooDecay::Flipped)");
  ws->factory("Decay::bkg4(Jpsi_Ct,lambdasym[0.69,0.02,5.0],resbkg,RooDecay::DoubleSided)");

  ws->factory("SUM::bkgPart1(fpm[0.95,0.,1.]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.,1.]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgctauTOT(fbkgTot[0.29,0.,1.]*resbkg,bkgPart2)");

  //psi' 

  ws->factory("Decay::bkg2P(Jpsi_Ct,lambdapP[0.42,0.05,1.5],resbkg,RooDecay::SingleSided)");
  ws->factory("Decay::bkg3P(Jpsi_Ct,lambdamP[0.79,0.02,1.5],resbkg,RooDecay::Flipped)");
  ws->factory("Decay::bkg4P(Jpsi_Ct,lambdasymP[0.69,0.02,5.0],resbkg,RooDecay::DoubleSided)");

  ws->factory("SUM::bkgPart1P(fpmP[0.95,0.,1.]*bkg2P,bkg3P)");
  ws->factory("SUM::bkgPart2P(fLivingP[0.9,0.,1.]*bkgPart1P,bkg4P)");
  ws->factory("SUM::bkgctauTOTP(fbkgTotP[0.29,0.,1.]*resbkg,bkgPart2P)");

  return;
}

void defineCTSignal(RooWorkspace *ws, 
		    RooDataHist *reducedNP, RooDataHist *reducedNPP)
{
  
  RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);  ws->import(sigNPW);
  RooHistPdfConv sigNPO("sigNPO","Non-prompt signal with outstanding gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigO")),*reducedNP);  ws->import(sigNPO);
  RooHistPdfConv sigNPM("sigNPM","Non-prompt signal with mastodontic gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigM")),*reducedNP);  ws->import(sigNPM);
  RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigN")),*reducedNP);   ws->import(sigNPN);

  // RooFormulaVar fracResNP("fracResNP","@0*@1",RooArgList(*(ws->var("fracRes")),*(ws->var("NpPrRatio"))));
  // ws->import(fracResNP);
  
  // RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPO")),*(ws->pdf("sigNPM")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->function("fracResNP")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));   ws->import(sigNP);
  RooAddPdf sigNP("sigNP","Non-prompt signal",RooArgSet(*(ws->pdf("sigNPW")),*(ws->pdf("sigNPO")),*(ws->pdf("sigNPM")),*(ws->pdf("sigNPN"))),RooArgSet(*(ws->var("fracRes")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));  ws->import(sigNP); 

  RooHistPdfConv sigNPWP("sigNPWP","Non-prompt signal with wide gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNPP);  ws->import(sigNPWP);
  RooHistPdfConv sigNPOP("sigNPOP","Non-prompt signal with outstanding gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigO")),*reducedNPP);  ws->import(sigNPOP);
  RooHistPdfConv sigNPMP("sigNPMP","Non-prompt signal with mastodontic gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigM")),*reducedNPP);  ws->import(sigNPMP);
  RooHistPdfConv sigNPNP("sigNPNP","Non-prompt signal with narrow gaussian",*(ws->var("Jpsi_Ct")),*(ws->var("meanResSigWP")),*(ws->var("sigmaResSigNP")),*reducedNPP);   ws->import(sigNPNP);

  // RooAddPdf sigNPP("sigNPP","Non-prompt signal 2",RooArgSet(*(ws->pdf("sigNPWP")),*(ws->pdf("sigNPOP")),*(ws->pdf("sigNPMP")),*(ws->pdf("sigNPNP"))),RooArgSet(*(ws->var("fracResP")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));  ws->import(sigNPP);
  RooAddPdf sigNPP("sigNPP","Non-prompt signal 2",RooArgSet(*(ws->pdf("sigNPWP")),*(ws->pdf("sigNPOP")),*(ws->pdf("sigNPMP")),*(ws->pdf("sigNPNP"))),RooArgSet(*(ws->var("fracRes")),*(ws->var("fracRes2")),*(ws->var("fracRes3"))));  ws->import(sigNPP);

  return;
}
*/

void getrange(string &varRange, float *varmin, float *varmax)
{
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
   cout << varRange.c_str() << ": range not valid!" << endl;
    assert(0);
  }

 return;
}

void setRanges(RooWorkspace *ws, float lmin, float lmax){

  const float JpsiMassMin = 2.5;
  const float JpsiMassMax = 3.5;

  // ws->var("Jpsi_CtTrue")->setRange(0.0,4.0);
  // ws->var("Jpsi_Ct")->setRange(-lmin,lmax);

  ws->var("Jpsi_Mass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("Jpsi_Mass")->setRange("left",JpsiMassMin,2.9);
  ws->var("Jpsi_Mass")->setRange("right",3.3,JpsiMassMax);

  // ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  // ws->cat("Jpsi_Type")->setRange("glbtrk","GT");

  // ws->cat("MCType")->setRange("prompt","PR");
  // ws->cat("MCType")->setRange("nonprompt","NP");
  // ws->cat("MCType")->setRange("bkg","BK");

  return;
}

/* RooBinning setMyBinning(float lmin, float lmax){

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
    rb2.addUniform(9,-lmin,-0.2);
    rb2.addUniform(40,-0.2,0.2);
    rb2.addUniform(16,0.2,0.5);
    rb2.addUniform(21,0.5,1.2);
    rb2.addUniform(8,1.2,lmax);
  }

  return rb2;
}
*/

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename, *filenameMC, *filenameMC2;
  int isGG = 2;
  string prange;
  string yrange;
  string lrange;
  bool prefitMass = false;
  float alphaF = 0.0;
  float enneF = 0.0;
  int syst = 0;

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

      case 's':
	syst = atoi(argv[i+1]);
	cout << "Systematics number " << syst << endl;
	break;	

	// case 'u':
	// prefitMass = true;
	// cout << "Will determine Bfrac, not NsigPR and NsigNP" << endl;
        // cout << "Signal MC mass pre-fitting activated" << endl;
	// break;

      case 'p':
	prange = argv[i+1];
	cout << "Range for pT is " << prange << " GeV/c" << endl;
        break;  
    
	// case 'l':
	// lrange = argv[i+1];
	// cout << "Range for l(J/psi) is -" << lrange << " mm" << endl;
        // break;  

      case 'y':
        yrange = argv[i+1];
        cout << "Range for |y| is " << yrange << endl;
        break;

      case 'a':
 	alphaF = atof(argv[i+1]);
        cout << "CB alpha is fixed to " << alphaF << endl;
 	break;

      case 'n':
 	enneF = atof(argv[i+1]);
        cout << "CB n is fixed to " << enneF << endl;
 	break;

//       case 'c':
//         filenameMC2 = argv[i+1];
//         cout << "File name for psiprime MC data is " << filenameMC2 << endl;
//         break;

//       case 'b':
// 	prefitBackground = true;
// 	cout << "The background ctau distribution will be prefitted and some parameters fixed" << endl;
// 	break; 
     
	// case 't': 
        //theTrigger = atoi(argv[i+1]);
        //cout << "Using trigger bit n. " << theTrigger << endl;
        //break; 
      } 

    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

//   TFile f2In(filenameMC);
//   f2In.cd();
//   RooDataSet *dataMC = (RooDataSet*)f2In.Get("dataJpsi");
//   dataMC->SetName("dataMC");

//   TFile f3In(filenameMC2);
//   f3In.cd();
//   RooDataSet *dataMC2 = (RooDataSet*)f3In.Get("dataPsip");
//   dataMC2->SetName("dataMC2");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("dataJpsi");
  data->SetName("data");

  float pmin, pmax; 
  float ymin, ymax;
  float lmin, lmax;

  // getrange(lrange,&lmin,&lmax);
  getrange(prange,&pmin,&pmax);
  getrange(yrange,&ymin,&ymax);
 
  // bool isTheSpecialBin = (fabs(ymin-1.6) < 0.001 && fabs(pmin-6.5) < 0.001);
  // if (isTheSpecialBin) cout << "Warning: for this bin we cheat in the figures to make reviewers happy" << endl;

  char reducestr[300];
  if (nVtx <= 0) {
    sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f", pmax,pmin,ymax,ymin);
    if (cutOnVert) {
      sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_Vprob > 0.01", pmax,pmin,ymax,ymin);
    }
    // TEST DELLA MINCHIA
    // sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_Ct > 1.5", pmax,pmin,ymax,ymin);  
  } else {
    sprintf(reducestr,"Jpsi_nPV == %d",nVtx);
    RooDataSet *tmpdata = (RooDataSet*)data->reduce(reducestr);
    cout << "All events with " << nVtx << " verteces = " << tmpdata->sumEntries() << endl;
    sprintf(reducestr,"Jpsi_Pt < %f && Jpsi_Pt > %f && abs(Jpsi_Y) < %f && abs(Jpsi_Y) > %f && Jpsi_nPV == %d", pmax,pmin,ymax,ymin,nVtx);
  }

//   RooDataSet *reddataMC = (RooDataSet*)dataMC->reduce(reducestr);
//   ws->import(*reddataMC);

//   RooDataSet *reddataMC2 = (RooDataSet*)dataMC2->reduce(reducestr);
//   ws->import(*reddataMC2);

  RooDataSet *reddata = (RooDataSet*)data->reduce(reducestr);
  ws->import(*reddata);

  setRanges(ws,lmin,lmax);

  string titlestr;
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");	
  ws->var("Jpsi_Mass")->SetTitle("#mu^{+} #mu^{-} invariant mass");
  // ws->var("Jpsi_Mass")->SetTitle("#psi(2S) mass");
  ws->var("Jpsi_Ct")->SetTitle("#font[12]{l}_{J/#psi}");
  
   // define binning for masses
  ws->var("Jpsi_Mass")->setBins(50);
  // ws->var("PsiP_Mass")->setBins(60);
  // ws->var("Jpsi_Ct")->setBins(45);

  // define binning for true lifetime
  // RooBinning rb(0.0001,4.0);
  // rb.addBoundary(-0.01);
  // rb.addUniform(100,0.0001,0.5);
  // rb.addUniform(15,0.5,1.0);
  // rb.addUniform(20,1.0,2.5);
  // rb.addUniform(5,2.5,4.0);
  // ws->var("Jpsi_CtTrue")->setBinning(rb);

  // define binning for lifetime
  // RooBinning rb2 = setMyBinning(lmin,lmax);
  // ws->var("Jpsi_Ct")->setBinning(rb2);

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

  RooDataHist *bindata = new RooDataHist("bindata","bindata",RooArgSet(*(ws->var("Jpsi_Mass"))),*reddata1);

//   RooDataSet *reddataJ = (RooDataSet*) reddata1->reduce("Jpsi_PsiP == Jpsi_PsiP::J");
//   RooDataSet *reddataP = (RooDataSet*) reddata1->reduce("Jpsi_PsiP == Jpsi_PsiP::P");
//   RooDataHist *bindataJ = new RooDataHist("bindataJ","bindataJ",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataJ);
//   RooDataHist *bindataP = new RooDataHist("bindataP","bindataP",RooArgSet(*(ws->var("PsiP_Mass")),*(ws->var("Jpsi_Ct"))),*reddataP);

  cout << "Number of events to fit  = " << bindata->sumEntries() << endl; 

  // Get subdatasets. Some of them are useful. Some, however, not  

  // RooDataSet *reddataSB = (RooDataSet*) reddata1->reduce("(Jpsi_PsiP == Jpsi_PsiP::J && (Jpsi_Mass < 2.9 || Jpsi_Mass > 3.3)) || (Jpsi_PsiP == Jpsi_PsiP::P && (PsiP_Mass < 3.45 || PsiP_Mass > 3.85))");
  RooDataSet *reddataSB = (RooDataSet*) reddata1->reduce("Jpsi_Mass < 2.9 || Jpsi_Mass > 3.3");

  // RooDataSet *reddataSBJ = (RooDataSet*) reddataSB->reduce("Jpsi_PsiP == Jpsi_PsiP::J");
  // RooDataSet *reddataSBP = (RooDataSet*) reddataSB->reduce("Jpsi_PsiP == Jpsi_PsiP::P");
  
  // cout << "Number of true events to fit  = " << reddataTr->sumEntries() << endl; 
  // RooDataHist* bindataPR = new RooDataHist("bindataPR","MC distribution for PR signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataPR);
  // RooDataHist* bindataPRP = new RooDataHist("bindataPRP","MC distribution for PR signal",RooArgSet(*(ws->var("PsiP_Mass")),*(ws->var("Jpsi_Ct"))),*reddataPRP);
  
  // RooDataHist* bindataNP = new RooDataHist("bindataNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataNP);

//   RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*reddataNP); 
//   RooDataHist* redMCNPP = new RooDataHist("redMCNPP","MC distribution for NP signal",RooArgSet(*(ws->var("Jpsi_CtTrue"))),*reddataNPP); 

//   // RooDataHist* bindataBK = new RooDataHist("bindataBK","MC distribution for background",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("Jpsi_Ct"))),*reddataBK);

//   RooDataHist* bindataSB = new RooDataHist("bindataSB","MC distribution for background",RooArgSet(*(ws->var("Jpsi_Mass")),*(ws->var("PsiP_Mass")),*(ws->var("Jpsi_Ct")),*(ws->cat("Jpsi_PsiP"))),*reddataSB);

  //JPSI MASS PARAMETRIZATION

  //background
  defineMassBackground(ws);

  //signal is the same for prompt and non-prompt
  defineMassSignal(ws);

  if (alphaF > 0.0001) {
    ws->var("alpha")->setVal(alphaF);
    ws->var("alpha")->setConstant(kTRUE);
  }
  if (enneF > 0.0001) {
    ws->var("enne")->setVal(enneF);
    ws->var("enne")->setConstant(kTRUE);
  }

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  // defineCTResol(ws);

  //background
  // defineCTBackground(ws);

  //signal
  // defineCTSignal(ws,redMCNP,redMCNPP);

  //putting all together
  // ws->factory("PROD::totBKG(expFunct,bkgctauTOT)");
  // ws->factory("PROD::totBKGP(expFunctP,bkgctauTOTP)");

  string partTit, partFile;
  if (isGG == 0) { partTit = " glb-glb "; partFile = "GG"; }
  else if (isGG == 1) { partTit = " glb-trk "; partFile = "GT"; }
  else { partTit = " all "; partFile = "ALL"; }

  int nFitPar;
  
  // ws->pdf("expFunct")->fitTo(*reddata1,Range("left,right"),SumW2Error(kTRUE));
  // ws->var("coefExp")->setConstant(kTRUE);
  // ws->pdf("sigCBGaussOneMean")->fitTo(*bindataTr,SumW2Error(kTRUE));
  // ws->var("enne")->setConstant(kTRUE);
  // ws->var("alpha")->setConstant(kTRUE);
  
  if (syst == 1) {
    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCBGaussOneMean,NBkg[2000.,10.,10000000.]*CPolFunct)");
  } else if (syst == 2) {
    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCB,NBkg[2000.,10.,10000000.]*expFunct)");
  } else {
    ws->factory("SUM::massPDF(NSig[5000.,10.,10000000.]*sigCBGaussOneMean,NBkg[2000.,10.,10000000.]*expFunct)");
  }

  RooFitResult* rfrmass; 
  if (pmax-pmin < 20.) {
    // Special bins
    if ( fabs(ymin-1.2) < 0.001 ) {
      if ( fabs(pmin-7.0) < 0.001 || fabs(pmin-7.5) < 0.001 || fabs(pmin-11.0) < 0.001 || fabs(pmin-12.0) < 0.001 || fabs(pmin-15.0) < 0.001 ) {
	ws->var("alpha")->setConstant(kTRUE);
	ws->var("enne")->setConstant(kTRUE);
	ws->var("sigmaSig1")->setMax(0.7);
      }
    }
    ws->pdf("massPDF")->fitTo(*reddata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2),Save(1));
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    rfrmass = ws->pdf("massPDF")->fitTo(*reddata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2),Save(1));
  } else {
    ws->pdf("massPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2),Save(1));
    ws->var("alpha")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    rfrmass = ws->pdf("massPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2),Save(1));
  }
  nFitPar = rfrmass->floatParsFinal().getSize() - 1;
  
  Double_t a = ws->var("NSig")->getVal();
  Double_t b = 0;
  const Double_t NSig_static[2] = {a,b}; 
   
  a = ws->var("NSig")->getError();  
  const Double_t Err_static[2] = {a,b};

  Double_t NBkg_static[2]; 
  NBkg_static[0] = ws->var("NBkg")->getVal();
  
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig2")->getVal();
  const double sigmaSig2 = ws->var("sigmaSig1")->getVal();
  double ecoeffGauss = ws->var("coeffGauss")->getError();
  if (ecoeffGauss > 0.3*coeffGauss) ecoeffGauss = 0.;
  const double esigmaSig1 = ws->var("sigmaSig2")->getError();
  const double esigmaSig2 = ws->var("sigmaSig1")->getError();

  float resol = sqrt(coeffGauss*sigmaSig1*sigmaSig1 + (1-coeffGauss)*sigmaSig2*sigmaSig2);
  float errresol = (0.5/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));	
  float bc = ws->var("coefExp")->getVal();
  float bSigReg =  ws->var("NBkg")->getVal()*(exp(3.146*0.02*bc)-exp(3.046*0.02*bc))/(exp(3.5*0.02*bc)-exp(2.6*0.02*bc));
  // float bSigReg = ws->var("NBkg")->getVal()/9.0;
  /*			 
  RooRealVar tempVar1("tempVar1","tempVar1",NSigNP_static[0]);
  RooRealVar tempVar2("tempVar2","tempVar2",NBkg_static[0]);
  RooRealVar tempVar3("tempVar3","tempVar3",NSigNP_static[1]);
  RooRealVar tempVar4("tempVar4","tempVar4",NBkg_static[1]);
   */
  // ### WRITE RESULTS
  cout << endl << "J/psi yields:" << endl;
  cout << "TOTAL Jpsi        : Fit : " << NSig_static[0] << " +/- " << Err_static[0] << endl;
  // cout << "PROMPT Jpsi       : Fit : " << NSigPR_static[0] << " +/- " << ErrPR_static[0] << endl;
  // cout << "NON-PROMPT Jpsi   : Fit : " << NSigNP_static[0] << " +/- " << ErrNP_static[0] << endl;
  // cout << "B fraction Jpsi   : Fit : " << Bfrac_static[0] << " +/- " << BfracErr_static[0] << endl;
  cout << endl << "psi(2S) yields:" << endl;
  cout << "TOTAL psi(2S)     : Fit : " << NSig_static[1] << " +/- " << Err_static[1] << endl;
  // cout << "PROMPT psi(2S)    : Fit : " << NSigPR_static[1] << " +/- " << ErrPR_static[1] << endl;
  // cout << "NON-PROMPT psi(2S): Fit : " << NSigNP_static[1] << " +/- " << ErrNP_static[1] << endl;
  // cout << "B fraction psi(2S): Fit : " << Bfrac_static[1] << " +/- " << BfracErr_static[1] << endl;
  // cout << endl << "Resolution Jpsi   : Fit : " << resol*1000. << " +/- " << errresol*1000. << " mum" << endl;
  // cout << "Resolution psi(2S): Fit : " << resolP*1000. << " +/- " << errresolP*1000. << " mum" << endl;
 
  char oFile[200];
  sprintf(oFile,"results/testVtxDen/results%s_pT%s_y%s.txt",partFile.c_str(),prange.c_str(),yrange.c_str());

  ofstream outputFile(oFile);
  outputFile << "TJ " << 0. << " " << NSig_static[0] << " " << Err_static[0] << endl;
  if (verbose) {
    outputFile << "RE " << 0. << " " << resol << " " << errresol << endl;
    outputFile << "SB " << 0. << " " << NSig_static[0]/bSigReg << endl;
  }
  // outputFile << "PJ " << 0. << " " << NSigPR_static[0] << " " << ErrPR_static[0] << endl;
  // outputFile << "NJ " << 0. << " " << NSigNP_static[0] << " " << ErrNP_static[0] << endl;
  // outputFile << "BJ " << 0. << " " << Bfrac_static[0] << " " << BfracErr_static[0] << endl;
  // outputFile << "TP " << 0. << " " << NSig_static[1] << " " << Err_static[1] << endl;
  // outputFile << "PP " << 0. << " " << NSigPR_static[1] << " " << ErrPR_static[1] << endl;
  // outputFile << "NP " << 0. << " " << NSigNP_static[1] << " " << ErrNP_static[1] << endl;
  // outputFile << "BP " << 0. << " " << Bfrac_static[1] << " " << BfracErr_static[1] << endl;
  outputFile << endl;


  /// ##### DRAW PLOTS ######
  // a) Jpsi mass / all mass

  RooPlot *mframe = ws->var("Jpsi_Mass")->frame();
  // mframe->SetMinimum(-500.); 
  RooHist *hresid, *hresidP;
  double chi2 = 0.;

  titlestr = " fit for" + partTit + "muons (J/ #psi mass projection), p_{T} = " + prange + " GeV/c and |y| = " + yrange;
  mframe->SetTitle(titlestr.c_str());

  reddata1->plotOn(mframe,DataError(RooAbsData::SumW2));

  if (syst == 2) {
    ws->pdf("massPDF")->plotOn(mframe,Components("CPolFunct"),LineStyle(kDashed),LineColor(kRed),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  } else {
    ws->pdf("massPDF")->plotOn(mframe,Components("expFunct"),LineStyle(kDashed),LineColor(kRed),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  }
  // RooAddPdf tempPDF("tempPDF","tempPDF",RooArgList(*(ws->pdf("sigCBGaussOneMean")),*(ws->pdf("expFunct"))),RooArgList(tempVar1,tempVar2));
  // tempPDF.plotOn(mframe,LineColor(kRed),Normalization(NSigNP_static[0] + NBkg_static[0],RooAbsReal::NumEvent));
  ws->pdf("massPDF")->plotOn(mframe,LineColor(kBlue),Normalization(reddata1->sumEntries(),RooAbsReal::NumEvent));
  hresid = mframe->pullHist();
  hresid->SetName("hresid");

  // NO RESIDUALS
  // TCanvas c1;
  // c1.cd();mframe->Draw();

  // WITH RESIDUALS
  TCanvas* c1 = new TCanvas("c1","The Canvas",200,10,440,600);
  c1->cd();

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();

  pad1->cd(); 
  // pad1->SetLogy(1); 
  mframe->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.045);
  // t->DrawLatex(0.7,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");
  t->DrawLatex(0.38,0.92,"CMS -  #sqrt{s} = 7 TeV"); 
  t->DrawLatex(0.38,0.86,"L = 36.7 pb^{-1}"); 

  Double_t fx[2], fy[2], fex[2], fey[2];
  TGraphErrors *gfake = new TGraphErrors(2,fx,fy,fex,fey);
  gfake->SetMarkerStyle(20);
  TH1F hfake1 = TH1F("hfake1","hfake1",100,200,300);
  hfake1.SetLineColor(kRed);
  hfake1.SetLineStyle(kDotted);
  hfake1.SetLineWidth(2);
  TH1F hfake2 = TH1F("hfake2","hfake2",100,200,300);
  hfake2.SetLineColor(kBlue);
  hfake2.SetLineWidth(2);

  TLegend * leg = new TLegend(0.22,0.68,0.58,0.835,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->AddEntry(gfake,"data","le1p");
  leg->AddEntry(&hfake2,"total fit","L");
  // if (superImpose) {
  //   leg->AddEntry(&hfake3,"bkgd + non-prompt","L");
  // } else {
  //  leg->AddEntry(&hfake4,"prompt","L");
  //  leg->AddEntry(&hfake3,"non-prompt","L");
  // }
  leg->AddEntry(&hfake1,"background","L");
  leg->Draw("same"); 

  // TLatex *t = new TLatex();
  // t->SetNDC();
  // t->SetTextAlign(22);
  // t->SetTextFont(63);
  // t->SetTextSize(0.05);
  // t->DrawLatex(0.6,0.9,"CMS Preliminary -  #sqrt{s} = 7 TeV");  

  RooPlot* mframeres =  ws->var("Jpsi_Mass")->frame(Title("Residuals Distribution")) ;
  mframeres->GetYaxis()->SetTitle("Pull");
  mframeres->SetLabelSize(0.08,"XYZ");
  mframeres->SetTitleSize(0.08,"XYZ");
  mframeres->SetTitleOffset(0.6,"Y");
  mframeres->SetTitleOffset(1.0,"X");
  mframeres->addPlotable(hresid,"P") ;  
  mframeres->SetMinimum(-4.0);
  mframeres->SetMaximum(-(mframeres->GetMinimum())); 

  pad2->cd(); mframeres->Draw();

  double *ypulls = hresid->GetY();
  unsigned int nBins = ws->var("Jpsi_Mass")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    cout << "Pull of bin " << i << " = " << ypulls[i] << endl;
    chi2 += ypulls[i]*ypulls[i]; 
    cout << "Partial chi2 = " << chi2 << endl;
    if (fabs(ypulls[i]) > 0.0001) nFullBins++;
  }
  // if (isTheSpecialBin) {
  //   hresid->SetPoint(1,-0.6,0.);
  //   hresid->SetPointError(1,0.,0.,0.,0.);
  // }
  cout << chi2 << endl;
  // chi2 /= (nFullBins - nFitPar);
  cout << chi2 << endl;
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
  } 

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  // t->SetTextFont(63);
  t2->SetTextSize(0.07);
  // sprintf(reducestr,"Reduced #chi^{2} = %f ; #chi^{2} probability = %f",chi2,TMath::Prob(chi2*nDOF,nDOF));
  cout << "Reduced chi2 = " << chi2 << endl;
  sprintf(reducestr,"#chi^{2}/n_{DoF} = %4.2f/%d",chi2,nFullBins - nFitPar);
  if (chi2 < 1000.) t2->DrawLatex(0.75,0.90,reducestr);
  
  c1->Update();

  // titlestr = "/afs/cern.ch/user/c/covarell/mynotes/tdr2/notes/AN-11-098/trunk/" + partFile + "massfitJpsi_pT" + prange + "_y" + yrange + ".pdf";
  titlestr = "pictures/" + partFile + "massfitJpsi_pT" + prange + "_y" + yrange + ".gif";
  c1->SaveAs(titlestr.c_str());

  return 1;
}
