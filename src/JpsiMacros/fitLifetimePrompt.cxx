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
  ws->factory("Gaussian::resGW(Jpsict,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.05,0.0001,0.5])");
  ws->factory("Gaussian::resGN(Jpsict,meanResSigN[0.,-0.1,0.1],sigmaResSigN[0.01,0.0001,0.5])");
  ws->factory("Gaussian::resGO(Jpsict,meanResSigW,sigmaResSigO[0.1,0.0001,0.3])");
  ws->factory("RooGExpModel::resLeft(Jpsict,sigmaResSigG[0.01,0.0001,0.5],resLife[0.1,0.0001,0.3])");
  // ws->factory("RooGExpModel::resRight(Jpsict,sigmaResSigG,resLife,RooGExpModel::Flipped)");

  // 3-Gaussian
  ws->factory("SUM::resol(fracRes[0.05,0.,1.]*resGW,fracRes2[0.05,0.,1.]*resGO,resGN)");
  // 2-Gaussian
  // ws->factory("SUM::resol(fracRes[0.05,0.,1.]*resGW,resGN)");
  // Gaussian + Gexp
  // ws->factory("SUM::resol(fracRes[0.5,0.,1.]*resLeft,resRight)");

  return;
}

void setRanges(RooWorkspace *ws){

  ws->var("Jpsict")->setRange(-0.5,0.5);

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
  if(isGG) c2.SaveAs("GGtimefit_Prompt_Lin.gif");
  else c2.SaveAs("GTtimefit_Prompt_Lin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit_Prompt_Log.gif");
  else c2.SaveAs("GTtimefit_Prompt_Log.gif");

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
