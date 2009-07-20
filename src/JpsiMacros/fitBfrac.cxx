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

  RooRealVar GGcoefside("c_GG_side","GG linear coefficient of bkg PDF",0.,-999.,999.);
  RooPolynomial GGsideFunct("GGsideFunc","GGsideFunc",JpsiMass,GGcoefside);

  RooRealVar meanSig1("meanSig1","Mean of the signal gaussian 1",3.1,2.8,3.2);
  RooRealVar sigmaSig1("sigmaSig1","#sigma of the signal gaussian 1",0.02,0.,1.);

  RooRealVar meanSig2("meanSig2","Mean of the signal gaussian 2",3.0,2.8,3.2);
  RooRealVar sigmaSig2("sigmaSig2","#sigma of the signal gaussian 2",0.04,0.,1.);

  RooGaussian signalG1("signalG1","Signal PDF 1",JpsiMass,meanSig1,sigmaSig1);
  RooGaussian signalG2("signalG2","Signal PDF 2",JpsiMass,meanSig2,sigmaSig2);

  RooRealVar coeffGauss("coeffGauss","Relative norm of the two signal gaussians",0.42,0.,1.);

  RooAddPdf sigPDF("sigPDF","Total signal pdf",signalG1,signalG2,coeffGauss);

  RooRealVar NSig("NSig","Number of signal events",5000.,10.,10000000.);
  RooRealVar NBkg("NBkg","Number of background events",1300.,10.,10000000.);

  RooAddPdf totPDF("totPDF","Total pdf",RooArgList(sigPDF,GGsideFunct),RooArgList(NSig,NBkg));

  totPDF.fitTo(*GGdata,Extended(1),Save(1));

  RooPlot *GGmframe = JpsiMass.frame();
  GGmframe->SetTitle("Mass fit for glb-glb muons");
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GGmframe);
  totPDF.plotOn(GGmframe,Components(RooArgSet(GGsideFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  totPDF.plotOn(GGmframe,Components(GGsideFunct),DrawOption("F"),FillColor(kRed));
  GGdata->plotOn(GGmframe,DataError(RooAbsData::SumW2));

  TCanvas c1;
  c1.cd();GGmframe->Draw();
  c1.SaveAs("GGmassfit.gif");

  //GT
  RooDataSet *GTdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT");
  GTdata->setWeightVar(MCweight);

  totPDF.fitTo(*GTdata,Extended(1),Save(1),Minos(0));

  RooPlot *GTmframe = JpsiMass.frame();
  GTmframe->SetTitle("Mass fit for glb-trk muons");
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GTmframe);
  totPDF.plotOn(GTmframe,Components(RooArgSet(GGsideFunct,sigPDF)),DrawOption("F"),FillColor(kGreen));
  totPDF.plotOn(GTmframe,Components(GGsideFunct),DrawOption("F"),FillColor(kRed));
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));

  TCanvas c2;
  c2.cd();GTmframe->Draw();
  c2.SaveAs("GTmassfit.gif");


  return 1;
}
