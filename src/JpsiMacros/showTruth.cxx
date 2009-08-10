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
#include "RooAddPdf.h"

using namespace RooFit;

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;

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
  const float JpsiMassMax = 3.5;
  const float JpsiCtMin = -1.0;
  const float JpsiCtMax = 5.0;

  RooRealVar JpsiMass("JpsiMass","#mu^{+}#mu^{-} mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  JpsiMass.setRange("left",2.6,2.9);
  JpsiMass.setRange("right",3.3,3.5);

  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,5.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  //GG
  RooDataSet *GGdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG");
  GGdata->setWeightVar(MCweight);
  RooDataSet *GGdataTr = (RooDataSet*)GGdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GGdataTr->setWeightVar(MCweight);
  RooDataSet *GGdataFa = (RooDataSet*)GGdata->reduce("MCType == MCType::BK");
  GGdataFa->setWeightVar(MCweight);

  RooPlot *GGmframe = JpsiMass.frame();
  GGmframe->SetTitle("Mass for glb-glb muons");
  GGdataTr->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  GGdataFa->plotOn(GGmframe,DataError(RooAbsData::SumW2),LineColor(kRed),MarkerColor(kRed));

  TCanvas c1;
  c1.cd();GGmframe->Draw();
  c1.SaveAs("GGmass.gif");

  //GT
  RooDataSet *GTdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT");
  GTdata->setWeightVar(MCweight);
  RooDataSet *GTdataTr = (RooDataSet*)GTdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GTdataTr->setWeightVar(MCweight);
  RooDataSet *GTdataFa = (RooDataSet*)GTdata->reduce("MCType == MCType::BK");
  GTdataFa->setWeightVar(MCweight);

  RooPlot *GTmframe = JpsiMass.frame();
  GTmframe->SetTitle("Mass for glb-trk muons");
  GTdataTr->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  GTdataFa->plotOn(GTmframe,DataError(RooAbsData::SumW2),LineColor(kRed),MarkerColor(kRed));

  TCanvas c2;
  c2.cd();GTmframe->Draw();
  c2.SaveAs("GTmass.gif");

  //GC
  RooDataSet *GCdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GC");
  GCdata->setWeightVar(MCweight);
  RooDataSet *GCdataTr = (RooDataSet*)GCdata->reduce("MCType == MCType::PR || MCType == MCType::NP");
  GCdataTr->setWeightVar(MCweight);
  RooDataSet *GCdataFa = (RooDataSet*)GCdata->reduce("MCType == MCType::BK");
  GCdataFa->setWeightVar(MCweight);

  RooPlot *GCmframe = JpsiMass.frame();
  GCmframe->SetTitle("Mass for glb-calo muons");
  GCdataTr->plotOn(GCmframe,DataError(RooAbsData::SumW2));
  GCdataFa->plotOn(GCmframe,DataError(RooAbsData::SumW2),LineColor(kRed),MarkerColor(kRed));

  TCanvas c3;
  c3.cd();GCmframe->Draw();
  c3.SaveAs("GCmass.gif");

  return 1;
}
