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
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;

void defineMassSignal(RooWorkspace *ws)
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

void defineMassBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(JpsiMass,{CoefPol1[-0.05,-1500.,1500.],CcoefPol2[0.]})");

  //Exponential
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-1.,-2.,0.1])");

  return;
}

void defineCTResol(RooWorkspace *ws)
{

  
  ws->factory("GaussModel::resGW(Jpsict,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.05,0.005,0.5])");
  ws->factory("GaussModel::resGN(Jpsict,meanResSigN[0.,-0.1,0.1],sigmaResSigN[0.01,0.005,0.5])");
  ws->factory("AddModel::resol({resGW,resGN},{fracRes[0.05,0.,1.]})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
  //bkg1 is the resolution function
  ws->factory("Decay::bkg2(Jpsict,lambdap[1.4,0.,5.],resol,RooDecay::SingleSided");
  ws->factory("Decay::bkg3(Jpsict,lambdam[1.88,0.,5.],resol,RooDecay::Flipped");
  ws->factory("Decay::bkg4(Jpsict,lambdasym[1.16,0.,10.],resol,RooDecay::DoubleSided");

  ws->factory("SUM::bkgPart1(fpm[1.,0.,1.]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.,1.]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgctauTOT(fbkgTot[0.5,0.,1.]*resol,bkgPart2)");

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

  RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);
  RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigN")),*(ws->var("sigmaResSigN")),*reducedNP);
  RooAddPdf sigNP("sigNP","Non-prompt signal",sigNPW,sigNPN,*(ws->var("fracRes")));
  // RooHistPdfConv sigNP("sigNP","Non-prompt signal with wide gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);

  ws->import(sigNP);

  return;
}

void setRanges(RooWorkspace *ws){

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;

  ws->var("JpsictTrue")->setRange(-1.0,5.0);
  // define binning
  RooBinning rb(-1.0,5.0);
  rb.addBoundary(-0.01);
  rb.addUniform(50,-0.01,0.5);
  rb.addUniform(10,0.5,1.0);
  rb.addUniform(15,1.0,3.0);
  rb.addUniform(2,3.0,5.0);
  ws->var("JpsictTrue")->setBinning(rb);

  ws->var("Jpsict")->setRange(-1.0,3.5);
  // define binning
  RooBinning rb2(-1.0,3.5);
  rb2.addBoundary(-0.5);
  rb2.addBoundary(-0.2);
  rb2.addBoundary(-0.1);
  rb2.addBoundary(-0.01);
  rb2.addUniform(22,-0.05,0.5);
  rb2.addUniform(10,0.5,1.0);
  rb2.addUniform(5,1.0,3.5);
  ws->var("Jpsict")->setBinning(rb2);

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

void drawResults(RooWorkspace *ws, const bool isGG)
{
  RooRealVar *JpsiMass = ws->var("JpsiMass");
  RooRealVar *Jpsict = ws->var("Jpsict");
  RooAbsPdf *totPDF = ws->pdf("totPDF");
  //RooAbsPdf *totsigPR = ws->pdf("totsigPR");
  //RooAbsPdf *totsigNP = ws->pdf("totsigNP");
  //RooAbsPdf *totBKG = ws->pdf("totBKG");

  RooPlot *mframe = JpsiMass->frame();

  if(isGG) mframe->SetTitle("2D fit for glb-glb muons (mass projection)");
  else mframe->SetTitle("2D fit for glb-trk muons (mass projection)");

  if(isGG) ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components("totsigNP,totBKG"),DrawOption("F"),FillColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components("totBKG"),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));

  if(isGG) ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c1;
  c1.cd();mframe->Draw();
  if(isGG) c1.SaveAs("2D_GGmassfit.gif");
  else c1.SaveAs("2D_GTmassfit.gif");

  RooPlot *tframe = Jpsict->frame();

  if(isGG) tframe->SetTitle("2D fit for glb-glb muons (c  #tau projection)");
  else tframe->SetTitle("2D fit for glb-trk muons (c  #tau projection)");

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(tframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,Components("totsigNP,totBKG"),DrawOption("F"),FillColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,Components("totBKG"),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(tframe,Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("2D_GGtimefit_Lin.gif");
  else c2.SaveAs("2DGTtimefit_Lin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("2D_GGtimefit_Log.gif");
  else c2.SaveAs("2D_GTtimefit_Log.gif");

  RooPlot *tframe1 = Jpsict->frame();

  if(isGG) tframe1->SetTitle("2D fit for glb-glb muons (c  #tau projection) - signal prompt");
  else tframe1->SetTitle("2D fit for glb-trk muons (c  #tau projection) - signal prompt");

  if(isGG) ws->data("data")->plotOn(tframe1,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::PR"));
  else ws->data("data")->plotOn(tframe1,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::PR"));

  totPDF->plotOn(tframe1,Components("totsigPR"),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c3;
  c3.cd();
  c3.cd();tframe1->Draw();
  if(isGG) c3.SaveAs("2D_GGTruePR_Lin.gif");
  else c3.SaveAs("2D_GTTruePR_Lin.gif");
  c3.SetLogy(1);
  c3.cd();tframe1->Draw();
  if(isGG) c3.SaveAs("2D_GGTruePR_Log.gif");
  else c3.SaveAs("2D_GTTruePR_Log.gif");

  RooPlot *tframe2 = Jpsict->frame();

  if(isGG) tframe2->SetTitle("2D fit for glb-glb muons (c  #tau projection) - signal non-prompt");
  else tframe2->SetTitle("2D fit for glb-trk muons (c  #tau projection) - signal non-prompt");

  if(isGG) ws->data("data")->plotOn(tframe2,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::NP"));
  else ws->data("data")->plotOn(tframe2,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::NP"));

  totPDF->plotOn(tframe2,Components("totsigNP"),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c4;
  c4.cd();
  c4.cd();tframe2->Draw();
  if(isGG) c4.SaveAs("2D_GGTrueNP_Lin.gif");
  else c4.SaveAs("2D_GTTrueNP_Lin.gif");
  c4.SetLogy(1);
  c4.cd();tframe2->Draw();
  if(isGG) c4.SaveAs("2D_GGTrueNP_Log.gif");
  else c4.SaveAs("2D_GTTrueNP_Log.gif");

  RooPlot *tframe3 = Jpsict->frame();

  if(isGG) tframe3->SetTitle("2D fit for glb-glb muons (c  #tau projection) - background");
  else tframe3->SetTitle("2D fit for glb-trk muons (c  #tau projection) - background");

  if(isGG) ws->data("data")->plotOn(tframe3,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"));
  else ws->data("data")->plotOn(tframe3,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"));

  totPDF->plotOn(tframe3,Components("totBKG"),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c5;
  c5.cd();
  c5.cd();tframe3->Draw();
  if(isGG) c5.SaveAs("2D_GGTrueBK_Lin.gif");
  else c5.SaveAs("2D_GTTrueBK_Lin.gif");
  c5.SetLogy(1);
  c5.cd();tframe3->Draw();
  if(isGG) c5.SaveAs("2D_GGTrueBK_Log.gif");
  else c5.SaveAs("2D_GTTrueBK_Log.gif");

  return;
}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  bool isGG;
  bool prefitSignalMass = false;
  bool prefitSignalCTau = false;
  bool prefitBackground = false;

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
	cout << "Is global global? flag is " << isGG << endl;
	break;

      case 'p':
	prefitSignalMass = true;
	cout << "The signal mass distribution will be prefitted on MC" << endl;
	break;

      case 'c':
	prefitSignalCTau = true;
	cout << "The signal ctau distribution will be prefitted on MC, values will be fixed in the background prefit (if any)" << endl;
	break;

      case 'b':
	prefitBackground = true;
	cout << "The background ctau distribution will be prefitted on MC (but not fixed!!!)" << endl;
	break;
      }
    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");

  data->setWeightVar("MCweight");

  ws->import(*data);

  setRanges(ws);

  ws->var("JpsiMass")->setBins(60);
  // ws->var("Jpsict")->setBins(45);

  //CONSIDER THE CASE
  RooDataSet *reddata1;

  if(isGG) reddata1 = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG");
  else reddata1 = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT");

  RooDataHist *bindata = new RooDataHist("bindata","bindata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict")),*(ws->cat("MCType"))),*reddata1);

  cout << "Number of events to fit  = " << bindata->sumEntries() << endl; 

  //get subdatasets. Some of them are useful. Some, however, not
  RooDataSet *reddataPR = (RooDataSet*) reddata1->reduce("MCType == MCType::PR");
  RooDataSet *reddataNP = (RooDataSet*) reddata1->reduce("MCType == MCType::NP");
  RooDataSet *reddataBK = (RooDataSet*) reddata1->reduce("MCType == MCType::BK");

  RooDataHist* bindataPR = new RooDataHist("bindataPR","MC distribution for PR signal",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataPR);

  RooDataHist* bindataNP = new RooDataHist("bindataNP","MC distribution for NP signal",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataNP);

  RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("JpsictTrue"))),*reddataNP); 

  RooDataHist* bindataBK = new RooDataHist("bindataBK","MC distribution for background",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict"))),*reddataBK);

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
  ws->factory("PROD::totsigPR(sigCBGauss,resol)");
  ws->factory("PROD::totsigNP(sigCBGauss,sigNP)");
  ws->factory("PROD::totBKG(expFunct,bkgctauTOT)");

  ws->factory("SUM::totPDF(NSigPR[4000.,10.,1000000.]*totsigPR,NSigNP[900.,10.,1000000.]*totsigNP,NBkg[1400.,10.,1000000.]*totBKG)");

  ws->loadSnapshot("fit2dpars_GG.root");

  if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*bindataPR,SumW2Error(kTRUE)/*,NumCPU(4)*/);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);
  }

  if(prefitSignalCTau){
    ws->pdf("resol")->fitTo(*bindataPR,SumW2Error(kTRUE)/*,NumCPU(4)*/);
  }

  if(prefitBackground){
    cout << "Prefitting background on " << bindataBK->sumEntries() << " MC events " << endl;

    if(prefitSignalCTau){
      //ws->var("meanResSigN")->setConstant(kTRUE);
      ws->var("meanResSigW")->setConstant(kTRUE);
      //ws->var("sigmaResSigN")->setConstant(kTRUE);
      ws->var("sigmaResSigW")->setConstant(kTRUE);
    }

    ws->pdf("bkgctauTOT")->fitTo(*bindataBK,SumW2Error(kTRUE)/*,NumCPU(4)*/);

    if(prefitSignalCTau){
      //ws->var("meanResSigN")->setConstant(kFALSE);
      ws->var("meanResSigW")->setConstant(kFALSE);
      //ws->var("sigmaResSigN")->setConstant(kFALSE);
      ws->var("sigmaResSigW")->setConstant(kFALSE);
    }

  }

  
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  ws->pdf("totPDF")->fitTo(*bindata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(4));

  ws->saveSnapshot("fit2dpart_GG.root",ws->components(),kFALSE);

  const Double_t NSigNP_static = ws->var("NSigNP")->getVal();
  const Double_t NSigPR_static = ws->var("NSigPR")->getVal();

  Double_t Bfrac = NSigNP_static/(NSigNP_static + NSigPR_static);
  cout << "B frac = " << Bfrac << " +/- " << Bfrac/NSigNP_static << endl;

  drawResults(ws,isGG);

  cout << endl << "J/psi yields:" << endl;
  cout << "PROMPT :     True MC : " << bindataPR->sumEntries() << " Fit : " << ws->var("NSigPR")->getVal() << " +/- " << ws->var("NSigPR")->getError() << endl;
  cout << "NON-PROMPT : True MC : " << bindataNP->sumEntries() << " Fit : " << ws->var("NSigNP")->getVal() << " +/- " << ws->var("NSigNP")->getError() << endl;

  return 1;
}
