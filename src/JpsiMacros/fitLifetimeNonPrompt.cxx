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
#include "RooUniformBinning.h"
#include "RooWorkspace.h"

#include "RooHistPdfConv.h"
#include "RooAddPdf.h"

using namespace RooFit;

//Convention: when necessary, use the following convention
// 0 Global-Global
// 1 Global-Tracker
// 2 Tracker-Tracker

void defineCTSignal(RooWorkspace *ws, const bool isGG)
{

  // VARIABLES
  RooRealVar *Jpsict = ws->var("Jpsict");
  RooRealVar *JpsictTrue = ws->var("JpsictTrue");
  RooRealVar meanResSigW("meanResSigW","Mean of the resolution wide gaussian",0.0004,-0.1,0.1);
  RooRealVar sigmaResSigW("sigmaResSigW","#sigma of the resolution wide gaussian",0.104,0.001,0.5);
  RooRealVar meanResSigN("meanResSigN","Mean of the resolution narrow gaussian",0.,-0.1,0.1);
  RooRealVar sigmaResSigN("sigmaResSigN","#sigma of the resolution narrow gaussian",0.034,0.001,0.5);
  RooRealVar fracRes("fracRes","Fraction of narrow/wider gaussians",0.172,0.,1.);
  
  // CT TRUE TEMPLATE
  RooDataSet* reddata2;
  if(isGG) reddata2 = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GG");
  else reddata2 = (RooDataSet*)ws->data("data")->reduce("JpsiType == JpsiType::GT");
  RooDataSet *reddataNP = (RooDataSet*)reddata2->reduce("MCType == MCType::NP");
  ws->var("JpsictTrue")->setRange(-1.0,5.0);
  // define binning
  RooBinning rb(-1.0,5.0);
  rb.addBoundary(-0.01);
  rb.addUniform(50,-0.01,0.5);
  rb.addUniform(30,0.5,1.0);
  rb.addUniform(15,1.0,3.0);
  rb.addUniform(2,3.0,5.0);
  ws->var("JpsictTrue")->setBinning(rb);

  RooDataHist* redMCNP = new RooDataHist("redMCNP","MC ctau distribution for NP signal",*(ws->var("JpsictTrue")),*reddataNP); 

  // RESOLUTION

  RooHistPdfConv histPdfW("histPdfW","histPdfW",*Jpsict,meanResSigW,sigmaResSigW,*redMCNP); 
  RooHistPdfConv histPdfN("histPdfN","histPdfN",*Jpsict,meanResSigN,sigmaResSigN,*redMCNP); 
  RooAddPdf histPdf("histPdf","histPdf",histPdfW,histPdfN,fracRes);

  ws->import(histPdf);

  // TEMPORARY PLOTS
  /* RooPlot *tframe = Jpsict->frame();
  int numEvents = 1000;
  ws->pdf("histPdf")->plotOn(tframe,Normalization(numEvents,RooAbsReal::NumEvent),LineWidth(1));
  ws->pdf("histPdf")->plotOn(tframe,Components("histPdfW"),LineColor(2),Normalization(numEvents,RooAbsReal::NumEvent),LineWidth(1));
  ws->pdf("histPdf")->plotOn(tframe,Components("histPdfN"),LineColor(1),Normalization(numEvents,RooAbsReal::NumEvent),LineWidth(1));
  
  TCanvas c1;
  c1.cd();
  tframe->Draw();
  c1.SaveAs("temp.gif");
  c1.SetLogy(1);
  c1.SaveAs("templog.gif");

  RooPlot *ttrueframe = JpsictTrue->frame();
  redMCNP->plotOn(ttrueframe);
  TCanvas c2;
  c2.cd();
  ttrueframe->Draw();
  
  c2.SaveAs("temp2.gif");
  c2.SetLogy(1);
  c2.SaveAs("temp2log.gif"); */

  RooRealVar JpsictRes("JpsictRes","J/psi ctau resolution",-0.4,0.4,"mm");

  const RooArgSet* thisRow = reddataNP->get();
  RooDataSet* dataRes = new RooDataSet("dataRes","Resolution",
                                        RooArgList(JpsictRes));

  for (Int_t iSamp = 0; iSamp < reddataNP->numEntries(); iSamp++)
    {
  
      thisRow = reddataNP->get(iSamp);

      RooRealVar* myJpsict = (RooRealVar*)thisRow->find("Jpsict");
      RooRealVar* myJpsictTrue = (RooRealVar*)thisRow->find("JpsictTrue");

      if ( myJpsictTrue->getVal() > 0.0001 ) {
	float Jpsires = myJpsict->getVal() - myJpsictTrue->getVal();
	JpsictRes.setVal(Jpsires);
	
	dataRes->add(RooArgSet(JpsictRes));
      }
    }

  RooPlot* frameRes = JpsictRes.frame();
  frameRes->SetTitle("c #tau resolution - MC non prompt");
  dataRes->plotOn(frameRes,DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2),MarkerSize(0.5)); 

  TCanvas c3;
  c3.cd();
  frameRes->Draw();
  
  c3.SaveAs("temp3.gif");
  c3.SetLogy(1);
  c3.SaveAs("temp3log.gif");

  return;
}

void setRanges(RooWorkspace *ws){

  // ws->var("Jpsict")->setRange(-1.0,1.0);

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

  if(isGG) tframe->SetTitle("Lifetime fit for non-prompt glb-glb J/   #psi");
  else tframe->SetTitle("Lifetime fit for non-prompt glb-trk J/   #psi");

  if(isGG) ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::NP && JpsiType == JpsiType::GG && JpsictTrue > 0.0001"));
  else ws->data("data")->plotOn(tframe,DataError(RooAbsData::SumW2),Cut("MCType == MCType::NP && JpsiType == JpsiType::GT && JpsictTrue > 0.0001"));

  ws->pdf("histPdf")->plotOn(tframe,Normalization(numEvents,RooAbsReal::NumEvent));

  TCanvas c2;
  c2.cd();
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit_nonPr_Lin.gif");
  else c2.SaveAs("GTtimefit_nonPr_Lin.gif");
  c2.SetLogy(1);
  c2.cd();tframe->Draw();
  if(isGG) c2.SaveAs("GGtimefit_nonPr_Log.gif");
  else c2.SaveAs("GTtimefit_nonPr_Log.gif");

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

  if(isGG) reddata1 = (RooDataSet*)data->reduce("MCType == MCType::NP && JpsiType == JpsiType::GG && JpsictTrue > 0.0001");
  else reddata1 = (RooDataSet*)data->reduce("MCType == MCType::NP && JpsiType == JpsiType::GT && JpsictTrue > 0.0001");

  RooDataHist *reddata = new RooDataHist("reddata","reddata",RooArgSet(*(ws->var("JpsiMass")),*(ws->var("Jpsict")),*(ws->cat("MCType"))),*reddata1);

  float numEvents = reddata1->sumEntries();
  cout << "Number of events to fit  = " << numEvents << endl; 

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  defineCTSignal(ws,isGG);

  ws->var("meanResSigN")->setConstant(kTRUE);
  ws->var("meanResSigW")->setConstant(kTRUE);

  //background
  ws->pdf("histPdf")->fitTo(*reddata,Minos(0),SumW2Error(kFALSE)/*,NumCPU(4)*/);
  // ws->pdf("histPdf")->fitTo(*reddata1,Minos(0),SumW2Error(kTRUE)/*,NumCPU(4)*/);

  drawResults(ws,isGG,numEvents); 

  return 1;
}
