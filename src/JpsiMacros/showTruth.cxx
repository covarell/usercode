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

void getrange(string &varRange, float *varmin, float *varmax){
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
    cout << varRange.c_str() << " not valid!" << endl;
    assert(0);
  }

 return;
}


int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  Int_t nevents = 0;
  string prange;
  string etarange;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){
      case 'f':{
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
 
      case 'p':{
	prange = argv[i+1];
	cout << "Range for pT is " << prange << " GeV/c" << endl;
        break;
      }
       
      case 'e':{
        etarange = argv[i+1];
        cout << "Range for |eta| is " << etarange << endl;
        break;
      }
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

  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
  RooRealVar* JpsictTrue = new RooRealVar("JpsictTrue","J/psi ctau true",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar JpsictRes("JpsictRes","J/psi ctau resolution",-0.4,0.4,"mm");
  JpsictRes.setBins(300);
  JpsictRes.setRange("fitRange",-0.25,0.25);
  JpsictRes.setRange("fitRange2",-0.35,0.35);

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,50.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  float pmin, pmax; 
  float etamin, etamax;

  getrange(prange,&pmin,&pmax);
  getrange(etarange,&etamin,&etamax);

  char reducestr[200];
  sprintf(reducestr,"abs(JpsictTrue) < 3.4999 && JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f", pmax,pmin,etamax,etamin);

  RooDataSet *OKdata = (RooDataSet*)data->reduce(RooArgSet(*Jpsict,*JpsictTrue,JpsiType,MCType),reducestr);
  const RooArgSet* thisRow = OKdata->get();

  RooDataSet* dataRes = new RooDataSet("dataRes","Resolution",
                                        RooArgList(JpsictRes));


  for (Int_t iSamp = 0; iSamp < OKdata->numEntries(); iSamp++)
    {
      if (nevents%100 == 0) cout << " >>> Processing event : " << nevents <<endl;
      nevents++;
      thisRow = OKdata->get(iSamp);

      RooRealVar* myJpsict = (RooRealVar*)thisRow->find("Jpsict");
      RooRealVar* myJpsictTrue = (RooRealVar*)thisRow->find("JpsictTrue");

      float Jpsires = myJpsict->getVal() - myJpsictTrue->getVal();
      JpsictRes.setVal(Jpsires);

      dataRes->add(RooArgSet(JpsictRes));

    }
  
  OKdata->merge(dataRes);
  cout << " Merging OK " <<endl;

  RooRealVar meanResSigW("meanResSigW","Mean of the resolution wide gaussian",-0.17,-1.,1.);
  RooRealVar sigmaResSigW("sigmaResSigW","#sigma of the resolution wide gaussian",0.02,0.001,5.);

  RooRealVar meanResSigN("meanResSigN","Mean of the resolution narrow gaussian",0.003,-1.,1.);
  RooRealVar sigmaResSigN("sigmaResSigN","#sigma of the resolution narrow gaussian",0.05,0.001,5.);

  RooRealVar meanResSigO("meanResSigO","Mean of the resolution outlier gaussian",-0.17,-1.,1.);
  RooRealVar sigmaResSigO("sigmaResSigO","#sigma of the resolution outlier gaussian",0.5,0.1,5.);

  RooGaussian resGW("resGW","Wide Gaussian resolution function",JpsictRes,meanResSigW,sigmaResSigW);
  RooGaussian resGN("resGN","Narrow Gaussian resolution function",JpsictRes,meanResSigN,sigmaResSigN);
  RooGaussian resGO("resGO","Outlier Gaussian resolution function",JpsictRes,meanResSigO,sigmaResSigO);

  RooRealVar fracRes("fracRes","Fraction of narrow/wider gaussians",0.13,0.,1.);
  RooRealVar fracRes2("fracRes2","Fraction of narrow/outlier gaussians",0.01,0.,1.);

  RooAddPdf resolGGPR("resolGGPR","resol",RooArgList(resGN,resGW,resGO),RooArgList(fracRes,fracRes2));
  RooAddPdf resolGGNP("resolGGNP","resol",RooArgList(resGN,resGW,resGO),RooArgList(fracRes,fracRes2));
  // RooAddPdf resol("resol","resol",RooArgList(resGW,resGN),RooArgList(fracRes));
  RooAddPdf resolGTPR("resolGTPR","resol",RooArgList(resGN,resGW,resGO),RooArgList(fracRes,fracRes2));
  RooAddPdf resolGTNP("resolGTNP","resol",RooArgList(resGN,resGW,resGO),RooArgList(fracRes,fracRes2));

  //GG
  TCanvas c1;
  c1.Divide(2,2);

  RooDataSet *mydataPRGG = (RooDataSet*)OKdata->reduce("JpsiType==JpsiType::GG && MCType==MCType::PR && abs(JpsictRes) < 0.399");
  RooDataHist histPRGG("histPRGG","GG prompt binned dataset",RooArgSet(JpsictRes),*mydataPRGG,1.0);
  resolGGPR.fitTo(histPRGG,Range("fitRange"));

  c1.cd(1);
  // gPad->SetLogy(1);
  RooPlot* frameGGPR = JpsictRes.frame();
  frameGGPR->SetTitle("c #tau resolution - MC prompt");
  histPRGG.plotOn(frameGGPR,DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2),MarkerSize(0.5)); 
  resolGGPR.plotOn(frameGGPR,Range("fitRange"));
  resolGGPR.plotOn(frameGGPR,Components(RooArgList(resGW,resGO)),LineStyle(kDashed),LineColor(4),Range("fitRange"));
  // resolGGPR.plotOn(frameGGPR,Components(resGO),LineStyle(kDashed),LineColor(4),Range("fitRange"));
  frameGGPR->Draw();

  RooDataSet *mydataNPGG = (RooDataSet*)OKdata->reduce("JpsiType==JpsiType::GG && MCType==MCType::NP && abs(JpsictRes) < 0.399");
  RooDataHist histNPGG("histNPGG","GG non-prompt binned dataset",RooArgSet(JpsictRes),*mydataNPGG,1.0);
  resolGGNP.fitTo(histNPGG,Range("fitRange"));

  c1.cd(2);
  // gPad->SetLogy(1);
  RooPlot* frameGGNP = JpsictRes.frame();
  frameGGNP->SetTitle("c #tau resolution - MC non-prompt");
  histNPGG.plotOn(frameGGNP,DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2),MarkerSize(0.5),RefreshNorm());
  resolGGNP.plotOn(frameGGNP,Range("fitRange"));
  resolGGNP.plotOn(frameGGNP,Components(RooArgList(resGW,resGO)),LineStyle(kDashed),LineColor(4),Range("fitRange"));
  // resolGGNP.plotOn(frameGGNP,Components(resGO),LineStyle(kDashed),LineColor(4),Range("fitRange"));
  frameGGNP->Draw();

  c1.cd(3);
  gPad->SetLogy(1);
  frameGGPR->Draw();
  /* RooPlot* frameGGBK = JpsictRes.frame();
  frameGGBK->SetTitle("c #tau resolution - MC background");
  OKdata->plotOn(frameGGBK,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),DataError(RooAbsData::SumW2));
  frameGGBK->Draw();*/

  c1.cd(4);
  gPad->SetLogy(1);
  frameGGNP->Draw();

  c1.SaveAs("GGctres.eps");

  //GT
  TCanvas c2;
  c2.Divide(2,2);

  RooDataSet *mydataPRGT = (RooDataSet*)OKdata->reduce("JpsiType==JpsiType::GT && MCType==MCType::PR && abs(JpsictRes) < 0.399");
  RooDataHist histPRGT("histPRGT","GT prompt binned dataset",RooArgSet(JpsictRes),*mydataPRGT,1.0);
  resolGTPR.fitTo(histPRGT,Range("fitRange2"));

  c2.cd(1);
  // gPad->SetLogy(1);
  RooPlot* frameGTPR = JpsictRes.frame();
  frameGTPR->SetTitle("c #tau resolution - MC prompt");
  histPRGT.plotOn(frameGTPR,DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2),MarkerSize(0.5)); 
  resolGTPR.plotOn(frameGTPR,Range("fitRange2"));
  resolGTPR.plotOn(frameGTPR,Components(RooArgList(resGW,resGO)),LineStyle(kDashed),LineColor(4),Range("fitRange2"));
  // resolGTPR.plotOn(frameGTPR,Components(resGO),LineStyle(kDashed),LineColor(4),Range("fitRange2"));
  frameGTPR->Draw();

  RooDataSet *mydataNPGT = (RooDataSet*)OKdata->reduce("JpsiType==JpsiType::GT && MCType==MCType::NP && abs(JpsictRes) < 0.399");
  RooDataHist histNPGT("histNPGT","GT non-prompt binned dataset",RooArgSet(JpsictRes),*mydataNPGT,1.0);
  resolGTNP.fitTo(histNPGT,Range("fitRange2"));

  c2.cd(2);
  // gPad->SetLogy(1);
  RooPlot* frameGTNP = JpsictRes.frame();
  frameGTNP->SetTitle("c #tau resolution - MC non-prompt");
  histNPGT.plotOn(frameGTNP,DataError(RooAbsData::SumW2),LineColor(2),MarkerColor(2),MarkerSize(0.5),RefreshNorm());
  resolGTNP.plotOn(frameGTNP,Range("fitRange2"));
  resolGTNP.plotOn(frameGTNP,Components(RooArgList(resGW,resGO)),LineStyle(kDashed),LineColor(4),Range("fitRange2"));
  // resolGTNP.plotOn(frameGTNP,Components(resGO),LineStyle(kDashed),LineColor(4),Range("fitRange2"));
  frameGTNP->Draw();

  c2.cd(3);
  gPad->SetLogy(1);
  frameGTPR->Draw();
  /* RooPlot* frameGTBK = JpsictRes.frame();
  frameGTBK->SetTitle("c #tau resolution - MC background");
  OKdata->plotOn(frameGTBK,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),DataError(RooAbsData::SumW2));
  frameGTBK->Draw();*/

  c2.cd(4);
  gPad->SetLogy(1);
  frameGTNP->Draw();

  c2.SaveAs("GTctres.eps");

  return 1;
}

