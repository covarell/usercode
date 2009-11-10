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
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooAddModel.h"
#include "RooGExpModel.h"
#include "RooFFTConvPdf.h"
#include "RooUniformBinning.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "RooExponential.h"

#include "RooHistPdfConv.h"

using namespace RooFit;

//list of parameters to save
//RooArgSet params;

void defineMassSignal(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // SIGNAL: 2-Gaussians
  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.09,3.05,3.15);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.03,0.,0.2);
  //params.add(meanSig1);
  //params.add(sigmaSig1);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.09,3.05,3.15);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.05,0.,0.2);
  //params.add(meanSig2);
  //params.add(sigmaSig2);

  // Different mean
  RooGaussian signalG1("signalG1","Signal PDF 1",*JpsiMass,meanSig1,sigmaSig1);
  RooGaussian signalG2("signalG2","Signal PDF 2",*JpsiMass,meanSig2,sigmaSig2);

  // Same mean
  RooGaussian signalG2OneMean("signalG2OneMean","Signal PDF 2",*JpsiMass,meanSig1,sigmaSig2);

  RooRealVar coeffGauss("coeffGauss","Relative norm of the two signal gaussians",0.8,0.,1.);
  //params.add(coeffGauss);

  // Different mean
  RooAddPdf sigPDF("sigPDF","Total signal pdf",signalG1,signalG2,coeffGauss);
  // Same mean
  RooAddPdf sigPDFOneMean("sigPDFOneMean","Total signal pdf",signalG1,signalG2OneMean,coeffGauss);

  // SIGNAL: Crystal Ball shape
  RooRealVar alpha("alpha","#alpha of CB",2.2,0.,5.);
  RooRealVar enne("enne","n of CB",0.33,0.,8.);
  //params.add(alpha);
  //params.add(enne);

  RooCBShape sigCB("sigCB","Signal CB PDF",*JpsiMass,meanSig1,sigmaSig1,alpha,enne);

  RooAddPdf sigCBGauss("sigCBGauss","Signal CB+Gauss PDF",sigCB,signalG2,coeffGauss);

  ws->import(RooArgSet(sigPDF,sigPDFOneMean,sigCBGauss),RecycleConflictNodes());

  return;
}

void defineMassBackground(RooWorkspace *ws){

  RooRealVar *JpsiMass = ws->var("JpsiMass");

  // BKG: first and second order polynomials
  RooRealVar CcoefPol1("CcoefPol1","linear coefficient of bkg PDF",-0.05,-1500.,1500.);
  RooRealVar CcoefPol2("CcoefPol2","quadratic coefficient of bkg PDF",0.1,-1.,1.);
  //params.add(CcoefPol1);
  //params.add(CcoefPol2);

  RooPolynomial CPolFunct1("CPolFunct","CPolFunct1",*JpsiMass,CcoefPol1);
  RooPolynomial CPolFunct2("CPolFunct2","CPolFunct2",*JpsiMass,RooArgList(CcoefPol1,CcoefPol2));
  
  CcoefPol2.setVal(0.0); 
  CcoefPol2.setConstant(kTRUE); 

  // BKG : exponential
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.);
  RooExponential expFunct("expFunct","expFunct",*JpsiMass,coefExp); 
  //params.add(coefExp);

  ws->import(RooArgSet(CPolFunct1,CPolFunct2,expFunct));

  return;
}

void defineCTResol(RooWorkspace *ws){

  RooRealVar *Jpsict = ws->var("Jpsict");

  RooRealVar meanResSigW("meanResSigW","Mean of the resolution wide gaussian",0.003,-1.,1.);
  RooRealVar sigmaResSigW("sigmaResSigW","#sigma of the resolution wide gaussian",0.05,0.001,5.);
  RooRealVar scaleK("scaleK","Scale factor of the resolution gaussian",1.,0.,10.);
  RooConstVar one("one","one",1.);
  //params.add(meanResSigW);
  //params.add(sigmaResSigW);
  //params.add(scaleK);

  scaleK.setConstant(kTRUE);

  RooRealVar meanResSigN("meanResSigN","Mean of the resolution narrow gaussian",-0.17,-1.,1.);
  RooRealVar sigmaResSigN("sigmaResSigN","#sigma of the resolution narrow gaussian",0.02,0.001,5.);
  //params.add(meanResSigN);
  //params.add(sigmaResSigN);

  RooGaussModel resGW("resGW","Wide Gaussian resolution function",*Jpsict,meanResSigW,sigmaResSigW,one,scaleK);
  RooGaussModel resGN("resGN","Narrow Gaussian resolution function",*Jpsict,meanResSigN,sigmaResSigN,one,scaleK);

  RooRealVar fracRes("fracRes","Fraction of narrow/wider gaussians",0.93,0.,1.);
  //params.add(fracRes);

  RooAddModel resol("resol","resol",RooArgList(resGW,resGN),RooArgList(fracRes));

  ws->import(resol);

  return;
}

void defineCTBackground(RooWorkspace *ws){

  RooRealVar *Jpsict = ws->var("Jpsict");

  RooRealVar lambdap("lambdap","tau of the positive background tail",0.4,0.,5.);
  RooRealVar lambdam("lambdam","tau of the negative background tail",3.88,0.,5.);
  RooRealVar lambdasym("lambdasym","tau of the symmetric background tail",0.16,0.,10.);
  //params.add(lambdap);
  //params.add(lambdam);
  //params.add(lambdasym);

  //bkg1 is the resolution function

  RooDecay bkg2("bkg2","One sided positive background",*Jpsict,lambdap,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::SingleSided);
  RooDecay bkg3("bkg3","One sided negative background",*Jpsict,lambdam,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::Flipped);
  RooDecay bkg4("bkg4","Symmetric background",*Jpsict,lambdasym,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::DoubleSided);

  //now, we could just compose the background together
  //but since we are not idiots, we compose them partially for stability
  RooRealVar fpm("fpm","Fraction of pos/neg tails",0.92,0.,1.);
  //params.add(fpm);

  RooAddPdf bkgPart1("bkgPart1","Sum of pos/neg backgrounds",bkg2,bkg3,fpm);

  RooRealVar fLiving("fLiving","Fraction of sym/asym living backgrounds",0.99,0.,1.);
  //params.add(fLiving);

  //FIXME
  RooAddPdf bkgPart2("bkgPart2","Sum of living backgrounds",bkgPart1,bkg4,fLiving);

  RooRealVar fbkgTot("fbkgTot","Fraction of delta living background",0.17,0.,1.);
  //params.add(fbkgTot);

  RooAddPdf bkgctauTOT("bkgctauTOT","Sum of all backgrounds",*(ws->pdf("resol")),bkgPart2,fbkgTot);

  ws->import(bkgctauTOT);

  return;
}

void defineCTSignal(RooWorkspace *ws, RooDataHist *reducedNP){

  //signal prompt, same as zero lifetime background

  //signal non-prompt
  //RooRealVar taueff("taueff","Effective tau of the B meson",0.348,0.05,1.);
  //params.add(taueff);

  //RooGExpModel physsigNP("physsigNP","Gauss + exp model",Jpsict,sigmaMC,taueff);

  //RooUniformBinning FFTbin(Jpsict.getMin(),Jpsict.getMax(),1000,"cache");
  //Jpsict.setBinning(FFTbin,"cache");
  //RooFFTConvPdf sigNP("sigNP","Non-prompt signal",Jpsict,physsigNP,resol);

  //RooDecay sigNP("sigNP","Non-prompt signal",*Jpsict,taueff,*(RooResolutionModel*)(ws->pdf("resol")),RooDecay::SingleSided);

  //get subdatasets. Some of them are useful. Some, however, not

  RooHistPdfConv sigNPW("sigNPW","Non-prompt signal with wide gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigW")),*(ws->var("sigmaResSigW")),*reducedNP);
  RooHistPdfConv sigNPN("sigNPN","Non-prompt signal with narrow gaussian",*(ws->var("Jpsict")),*(ws->var("meanResSigN")),*(ws->var("sigmaResSigN")),*reducedNP);

  RooAddPdf sigNP("sigNP","Non-prompt signal",sigNPW,sigNPN,*(ws->var("fracRes")));

  ws->import(sigNP);

  return;
}

void setRanges(RooWorkspace *ws){

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;
  //const float JpsiCtMin = -1.0;
  //const float JpsiCtMax = 5.0;

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

void drawResults(RooWorkspace *ws, const bool isGG){

  RooRealVar *JpsiMass = ws->var("JpsiMass");
  RooRealVar *Jpsict = ws->var("Jpsict");
  RooAbsPdf *totPDF = ws->pdf("totPDF");
  RooAbsPdf *totsigPR = ws->pdf("totsigPR");
  RooAbsPdf *totsigNP = ws->pdf("totsigNP");
  RooAbsPdf *totBKG = ws->pdf("totBKG");


  RooPlot *mframe = JpsiMass->frame();

  if(isGG) mframe->SetTitle("2D fit for glb-glb muons (mass projection)");
  else mframe->SetTitle("2D fit for glb-trk muons (mass projection)");

  if(isGG) ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components(RooArgSet(*totsigNP,*totBKG)),DrawOption("F"),FillColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(mframe,Components(RooArgSet(*totBKG)),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));

  if(isGG) ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(mframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(mframe,Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c1;
  c1.cd();mframe->Draw();
  if(isGG) c1.SaveAs("2DGGmassfit.gif");
  else c1.SaveAs("2DGTmassfit.gif");

  RooPlot *tframe = Jpsict->frame();

  if(isGG) tframe->SetTitle("2D fit for glb-glb muons (c  #tau projection)");
  else tframe->SetTitle("2D fit for glb-trk muons (c  #tau projection)");

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(tframe,Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,DrawOption("F"),FillColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,Components(RooArgSet(*totsigNP,*totBKG)),DrawOption("F"),FillColor(kBlue),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->plotOn(tframe,Components(RooArgSet(*totBKG)),DrawOption("F"),FillColor(kRed),Normalization(1.0,RooAbsReal::RelativeExpected));

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT"));

  totPDF->plotOn(tframe,Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("2DGGtimefitLin.gif");
  else c2.SaveAs("2DGTtimefitLin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("2DGGtimefit.gif");
  else c2.SaveAs("2DGTtimefit.gif");

  RooPlot *tframe1 = Jpsict->frame();

  if(isGG) tframe1->SetTitle("2D fit for glb-glb muons (c  #tau projection) - signal prompt");
  else tframe1->SetTitle("2D fit for glb-trk muons (c  #tau projection) - signal prompt");

  if(isGG) ws->data("data")->plotOn(tframe1,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::PR"));
  else ws->data("data")->plotOn(tframe1,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::PR"));

  totPDF->plotOn(tframe1,Components(RooArgSet(*totsigPR)),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c3;
  c3.cd();
  c3.cd();tframe1->Draw();
  if(isGG) c3.SaveAs("2DGGTruePRLin.gif");
  else c3.SaveAs("2DGTTruePRLin.gif");
  c3.SetLogy(1);
  c3.cd();tframe1->Draw();
  if(isGG) c3.SaveAs("2DGGTruePR.gif");
  else c3.SaveAs("2DGTTruePR.gif");

  RooPlot *tframe2 = Jpsict->frame();

  if(isGG) tframe2->SetTitle("2D fit for glb-glb muons (c  #tau projection) - signal non-prompt");
  else tframe2->SetTitle("2D fit for glb-trk muons (c  #tau projection) - signal non-prompt");

  if(isGG) ws->data("data")->plotOn(tframe2,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::NP"));
  else ws->data("data")->plotOn(tframe2,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::NP"));

  totPDF->plotOn(tframe2,Components(RooArgSet(*totsigNP)),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c4;
  c4.cd();
  c4.cd();tframe2->Draw();
  if(isGG) c4.SaveAs("2DGGTrueNPLin.gif");
  else c4.SaveAs("2DGTTrueNPLin.gif");
  c4.SetLogy(1);
  c4.cd();tframe2->Draw();
  if(isGG) c4.SaveAs("2DGGTrueNP.gif");
  else c4.SaveAs("2DGTTrueNP.gif");

  RooPlot *tframe3 = Jpsict->frame();

  if(isGG) tframe3->SetTitle("2D fit for glb-glb muons (c  #tau projection) - background");
  else tframe3->SetTitle("2D fit for glb-trk muons (c  #tau projection) - background");

  if(isGG) ws->data("data")->plotOn(tframe3,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"));
  else ws->data("data")->plotOn(tframe3,DataError(RooAbsData::SumW2),Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"));

  totPDF->plotOn(tframe3,Components(RooArgSet(*totBKG)),LineColor(kGreen),Normalization(1.0,RooAbsReal::RelativeExpected));

  TCanvas c5;
  c5.cd();
  c5.cd();tframe3->Draw();
  if(isGG) c5.SaveAs("2DGGTrueBKLin.gif");
  else c5.SaveAs("2DGTTrueBKLin.gif");
  c5.SetLogy(1);
  c5.cd();tframe3->Draw();
  if(isGG) c5.SaveAs("2DGGTrueBK.gif");
  else c5.SaveAs("2DGTTrueBK.gif");

  return;
}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  bool isGG;
  bool prefitSignalMass = false;
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

      case 'b':
	prefitBackground = true;
	cout << "The background ct distribution will be prefitted on MC (but not fixed!!!)" << endl;
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

  ws->var("JpsiMass")->setBins(25);
  ws->var("Jpsict")->setBins(25);

  //CONSIDER THE CASE
  RooDataSet *reddata;

  if(isGG) reddata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG");
  else reddata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT");

  //RooDataHist *reddata = new RooDataHist("reddata","reddata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict")),*(ws->cat("MCType"))),*reddata1);

  cout << "Number of events to fit  = " << reddata->numEntries(kTRUE) << endl; 

  //get subdatasets. Some of them are useful. Some, however, not
  RooDataSet *reddataPR = (RooDataSet*) reddata->reduce("MCType == MCType::PR");
  RooDataSet *reddataNP = (RooDataSet*) reddata->reduce("MCType == MCType::NP");
  RooDataSet *reddataBK = (RooDataSet*) reddata->reduce("MCType == MCType::BK");

  RooDataHist* redMCNP = new RooDataHist("redMCNP","MC distribution for NP signal",RooArgSet(*(ws->var("JpsictTrue"))),*reddataNP); 

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
  RooProdPdf totsigPR("totsigPR","Total prompt signal",RooArgList(*(ws->pdf("sigCBGauss")),*(ws->pdf("resol"))));
  RooProdPdf totsigNP("totsigNP","Total non-prompt signal",RooArgList(RooArgList(*(ws->pdf("sigCBGauss")),*(ws->pdf("sigNP")))));
  RooProdPdf totBKG("totBKG","Total background",RooArgList(*(ws->pdf("expFunct")),*(ws->pdf("bkgctauTOT"))));

  RooRealVar NSigPR("NSigPR","Number of prompt signal events",4000.,10.,1000000.);
  RooRealVar NSigNP("NSigNP","Number of non-prompt signal events",900.,10.,1000000.);
  RooRealVar NBkg("NBkg","Number of background events",1400.,10.,1000000.);
  //params.add(NSigPR);
  //params.add(NSigNP);
  //params.add(NBkg);

  RooAddPdf totPDF("totPDF","Total PDF",RooArgList(totsigPR,totsigNP,totBKG),RooArgList(NSigPR,NSigNP,NBkg));

  ws->import(totPDF);

  ws->loadSnapshot("fit2dpars_GG.root");

  if(prefitSignalMass){
    ws->pdf("sigCBGauss")->fitTo(*reddataPR,NumCPU(4),Minos(1));
    ws->var("enne")->setConstant(kTRUE);
    ws->var("alpha")->setConstant(kTRUE);
  }

  if(prefitBackground){
    cout << "Prefitting background on " << reddataBK->numEntries(kTRUE) << " MC events " << endl;
    ws->pdf("bkgctauTOT")->fitTo(*reddataBK/*,NumCPU(4)*/);
    //ws->pdf("resol")->fitTo(*reddataBK/*,NumCPU(4)*/);
  }

  ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Save(1),Minos(0),NumCPU(4));

  ws->saveSnapshot("fit2dpart_GG.root",ws->components(),kFALSE);

  Double_t Bfrac = NSigNP.getVal()/(NSigNP.getVal() + NSigPR.getVal());
  cout << "B frac = " << Bfrac << " +/- " << Bfrac/NSigNP.getVal() << endl;

  drawResults(ws,isGG);

  cout << endl << "GG J/psi yields:" << endl;
  cout << "PROMPT :     True MC : " << reddataPR->numEntries(true) << " Fit : " << ws->var("NSigPR")->getVal() << " +/- " << ws->var("NSigPR")->getError() << endl;
  cout << "NON-PROMPT : True MC : " << reddataNP->numEntries(true) << " Fit : " << ws->var("NSigNP")->getVal() << " +/- " << ws->var("NSigNP")->getError() << endl;

  return 1;
}
