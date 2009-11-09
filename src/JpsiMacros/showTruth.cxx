// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"

using namespace RooFit;

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  Int_t nevents = 0;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){
      case 'f':
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
    }
    }
  }

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");

  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.6;
  const float JpsiCtMin = -1.0;
  const float JpsiCtMax = 3.5;

  RooRealVar JpsictRes("JpsictRes","J/psi ctau resolution",-2.0,2.0,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,50.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  RooDataSet *OKdata = (RooDataSet*)data->reduce("JpsictTrue > -10.0");
  const RooArgSet* thisRow = OKdata->get();

  RooDataSet* dataRes = new RooDataSet("dataRes","Resolution",
                                        RooArgList(JpsictRes));


  for (Int_t iSamp = 0; iSamp < OKdata->numEntries(); iSamp++)
    {
      if (nevents%100 == 0) cout << " >>> Processing event : " << nevents <<endl;
      nevents++;
      thisRow = OKdata->get(iSamp);

      RooRealVar* Jpsict = (RooRealVar*)thisRow->find("Jpsict");
      RooRealVar* JpsictTrue = (RooRealVar*)thisRow->find("JpsictTrue");

      float Jpsires = Jpsict->getVal() - JpsictTrue->getVal();
      JpsictRes.setVal(Jpsires);

      dataRes->add(RooArgSet(JpsictRes));

    }
  
  OKdata->merge(dataRes);
  cout << " Merging OK " <<endl;

  //GG
  TCanvas c1;
  c1.Divide(2,2);

  c1.cd(1);
  gPad->SetLogy(1);
  RooPlot* frameGGPR = JpsictRes.frame();
  frameGGPR->SetTitle("c #tau resolution - MC prompt");
  OKdata->plotOn(frameGGPR,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::PR"),DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2));
  frameGGPR->Draw();

  c1.cd(2);
  gPad->SetLogy(1);
  RooPlot* frameGGNP = JpsictRes.frame();
  frameGGNP->SetTitle("c #tau resolution - MC non-prompt");
  OKdata->plotOn(frameGGNP,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::NP"),DataError(RooAbsData::SumW2),LineColor(4),MarkerColor(4));
  frameGGNP->Draw();

  c1.cd(3);
  gPad->SetLogy(1);
  RooPlot* frameGGBK = JpsictRes.frame();
  frameGGBK->SetTitle("c #tau resolution - MC background");
  OKdata->plotOn(frameGGBK,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),DataError(RooAbsData::SumW2));
  frameGGBK->Draw();

  c1.SaveAs("GGctres.eps");

  //GT
  TCanvas c2;
  c2.Divide(2,2);

  c2.cd(1);
  gPad->SetLogy(1);
  RooPlot* frameGTPR = JpsictRes.frame();
  frameGTPR->SetTitle("c #tau resolution - MC prompt");
  OKdata->plotOn(frameGTPR,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::PR"),DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2));
  frameGTPR->Draw();

  c2.cd(2);
  gPad->SetLogy(1);
  RooPlot* frameGTNP = JpsictRes.frame();
  frameGTNP->SetTitle("c #tau resolution - MC non-prompt");
  OKdata->plotOn(frameGTNP,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::NP"),DataError(RooAbsData::SumW2),LineColor(4),MarkerColor(4));
  frameGTNP->Draw();

  c2.cd(3);
  gPad->SetLogy(1);
  RooPlot* frameGTBK = JpsictRes.frame();
  frameGTBK->SetTitle("c #tau resolution - MC background");
  OKdata->plotOn(frameGTBK,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),DataError(RooAbsData::SumW2));
  frameGTBK->Draw();

  c2.SaveAs("GTctres.eps");

  return 1;
}
