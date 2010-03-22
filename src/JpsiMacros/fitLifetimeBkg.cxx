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
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"

using namespace RooFit;

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

void defineCTResol(RooWorkspace *ws)
{

  ws->factory("GaussModel::resGW(Jpsict,meanResSigW[0.,-0.1,0.1],sigmaResSigW[0.05,0.005,0.5])");
  ws->factory("GaussModel::resGN(Jpsict,meanResSigN[0.,-0.1,0.1],sigmaResSigN[0.01,0.005,0.5])");
  ws->factory("AddModel::resol({resGW,resGN},{fracRes[0.05,0.,1.]})");

  return;
}

void defineCTBackground(RooWorkspace *ws)
{
  ws->factory("Decay::bkg2(Jpsict,lambdap[1.4,0.,5.],resol,RooDecay::SingleSided");
  ws->factory("Decay::bkg3(Jpsict,lambdam[1.88,0.,5.],resol,RooDecay::Flipped");
  ws->factory("Decay::bkg4(Jpsict,lambdasym[1.16,0.,10.],resol,RooDecay::DoubleSided");

  ws->factory("SUM::bkgPart1(fpm[1.,0.,1.]*bkg2,bkg3)");
  ws->factory("SUM::bkgPart2(fLiving[0.9,0.,1.]*bkgPart1,bkg4)");
  ws->factory("SUM::bkgctauTOT(fbkgTot[0.5,0.,1.]*resol,bkgPart2)");

  return;
}

void setRanges(RooWorkspace *ws){

  ws->var("Jpsict")->setRange(-1.0,3.5);

  ws->cat("JpsiType")->setRange("glbglb","GG");
  ws->cat("JpsiType")->setRange("glbtrk","GT");

  ws->cat("MCType")->setRange("prompt","PR");
  ws->cat("MCType")->setRange("nonprompt","NP");
  ws->cat("MCType")->setRange("bkg","BK");

  return;
}

void drawResults(RooWorkspace *ws, const bool isGG, float numEvents, RooBinning binning)
{

  RooRealVar *Jpsict = ws->var("Jpsict");

  RooPlot *tframe = Jpsict->frame();

  if(isGG) tframe->SetTitle("Lifetime fit for background glb-glb J/  #psi");
  else tframe->SetTitle("Lifetime fit for background glb-trk J/  #psi");

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(binning),Cut("MCType == MCType::BK && JpsiType == JpsiType::GG"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(binning),Cut("MCType == MCType::BK && JpsiType == JpsiType::GT"));

  ws->pdf("bkgctauTOT")->plotOn(tframe,Normalization(numEvents,RooAbsReal::NumEvent));

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit_Bkg_Lin.gif");
  else c2.SaveAs("GTtimefit_Bkg_Lin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit_Bkg_Log.gif");
  else c2.SaveAs("GTtimefit_Bkg_Log.gif");

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

  // ws->var("Jpsict")->setBins(80);
  RooBinning rb2(-1.0,3.5);
  rb2.addBoundary(-0.5);
  rb2.addBoundary(-0.2);
  rb2.addUniform(30,-0.1,0.5);
  rb2.addUniform(10,0.5,1.0);
  rb2.addUniform(5,1.0,3.5);
  ws->var("Jpsict")->setBinning(rb2);

  //CONSIDER THE CASE
  RooDataSet *reddata1;

  if(isGG) reddata1 = (RooDataSet*)data->reduce("MCType == MCType::BK && JpsiType == JpsiType::GG");
  else reddata1 = (RooDataSet*)data->reduce("MCType == MCType::BK && JpsiType == JpsiType::GT");

  RooDataHist *reddata = new RooDataHist("reddata","reddata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict")),*(ws->cat("MCType"))),*reddata1);

  float numEvents = reddata1->sumEntries();
  cout << "Number of events to fit  = " << numEvents << endl; 

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  defineCTResol(ws);

  //background
  defineCTBackground(ws);

  ws->var("meanResSigN")->setConstant(kTRUE);
  ws->var("meanResSigW")->setConstant(kTRUE);
  ws->var("fpm")->setConstant(kTRUE);
  ws->var("fLiving")->setConstant(kTRUE);

  ws->pdf("bkgctauTOT")->fitTo(*reddata,Minos(0),SumW2Error(kFALSE)/*,NumCPU(4)*/);
  // ws->pdf("bkgctauTOT")->fitTo(*reddata1,Minos(0),SumW2Error(kTRUE)/*,NumCPU(4)*/);

  //  ws->saveSnapshot("fit2dpart_GG.root",ws->components(),kFALSE);

  drawResults(ws,isGG,numEvents,rb2);

  return 1;
}
