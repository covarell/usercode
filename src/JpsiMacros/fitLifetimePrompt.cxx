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
//#include "RooUniformBinning.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

void defineCTSignal(RooWorkspace *ws)
{
  ws->factory("Gaussian::resGW(Jpsict,meanResSigW[0.003,-1.,1.],sigmaResSigW[0.05,0.001,5.])");
  ws->factory("Gaussian::resGN(Jpsict,meanResSigN[0.003,-1.,1.],sigmaResSigN[0.05,0.001,5.])");
  ws->factory("SUM::resol(fracRes[0.5,0.,1.]*resGW,resGN)");

  return;
}

void setRanges(RooWorkspace *ws){

  ws->var("Jpsict")->setRange(-1.0,1.0);

  ws->cat("JpsiType")->setRange("glbglb","GG");
  ws->cat("JpsiType")->setRange("glbtrk","GT");

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

void drawResults(RooWorkspace *ws, const bool isGG, float numEvents)
{

  RooRealVar *Jpsict = ws->var("Jpsict");

  RooPlot *tframe = Jpsict->frame();

  if(isGG) tframe->SetTitle("Lifetime fit for prompt glb-glb J/  #psi");
  else tframe->SetTitle("Lifetime fit for prompt glb-trk J/  #psi");

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::PR && JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::PR && JpsiType == JpsiType::GT"));

  ws->pdf("resol")->plotOn(tframe,Normalization(numEvents,RooAbsReal::NumEvent));

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefitLin.gif");
  else c2.SaveAs("GTtimefitLin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit.gif");
  else c2.SaveAs("GTtimefit.gif");

  return;
}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  bool isGG;

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

  // ws->var("JpsiMass")->setBins(35);
  ws->var("Jpsict")->setBins(50);

  //CONSIDER THE CASE
  RooDataSet *reddata1;

  if(isGG) reddata1 = (RooDataSet*)data->reduce("MCType == MCType::PR && JpsiType == JpsiType::GG");
  else reddata1 = (RooDataSet*)data->reduce("MCType == MCType::PR && JpsiType == JpsiType::GT");

  RooDataHist *reddata = new RooDataHist("reddata","reddata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict")),*(ws->cat("MCType"))),*reddata1);

  float numEvents = reddata1->sumEntries();
  cout << "Number of events to fit  = " << numEvents << endl; 

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  defineCTSignal(ws);

  //background
  ws->pdf("resol")->fitTo(*reddata,Minos(0),SumW2Error(kFALSE)/*,NumCPU(4)*/);
  // ws->pdf("resol")->fitTo(*reddata1,Minos(0),SumW2Error(kTRUE)/*,NumCPU(4)*/);

  //  ws->saveSnapshot("fit2dpart_GG.root",ws->components(),kFALSE);

  drawResults(ws,isGG,numEvents);

  return 1;
}
