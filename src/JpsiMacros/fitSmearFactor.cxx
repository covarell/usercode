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

  RooRealVar Jpsict("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,5.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  //CONSIDER THE GG CASE
  RooDataSet *GGdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG && MCType == MCType::NP");
  GGdata->setWeightVar(MCweight);

  //JPSI CTAU PARAMETRIZATION

  RooRealVar meanResSigW("meanResSigW","Mean of the resolution wide gaussian",-0.000961866,-1.,1.);
  RooRealVar sigmaResSigW("sigmaResSigW","#sigma of the resolution wide gaussian",0.029023,0.,5.);
  RooRealVar scaleK("scaleK","Scale factor of the resolution gaussian",1.,0.,10.);
  RooConstVar one("one","one",1.);

  scaleK.setConstant(kTRUE);
  meanResSigW.setConstant(kTRUE);
  sigmaResSigW.setConstant(kTRUE);

  RooRealVar meanResSigN("meanResSigN","Mean of the resolution narrow gaussian",0.0027155,-1.,1.);
  RooRealVar sigmaResSigN("sigmaResSigN","#sigma of the resolution narrow gaussian",0.065901,0.,5.);

  meanResSigN.setConstant(kTRUE);
  sigmaResSigN.setConstant(kTRUE);

  RooGaussModel resGW("resGW","Wide Gaussian resolution function",Jpsict,meanResSigW,sigmaResSigW,one,scaleK);
  RooGaussModel resGN("resGN","Narrow Gaussian resolution function",Jpsict,meanResSigN,sigmaResSigN,one,scaleK);

  RooRealVar fracRes("fracRes","Fraction of narrow/wider gaussians",0.67362,0.,1.);

  fracRes.setConstant(kTRUE);

  RooAddModel resol("resol","resol",RooArgList(resGW,resGN),RooArgList(fracRes));

  RooRealVar taueff("taueff","Effective tau of the B meson",0.3,0.,1.);

  RooDecay sigNP("sigNP","Non-prompt signal",Jpsict,taueff,resol,RooDecay::SingleSided);

  RooFitResult* fitRes = sigNP.fitTo(*GGdata,Save(1),Minos(0));

  RooPlot *GGtframe = Jpsict.frame();
  GGdata->plotOn(GGtframe,DataError(RooAbsData::SumW2));
  sigNP.plotOn(GGtframe);

  cout << "chi2 = " << GGtframe->chiSquare();

  TCanvas c2;
  c2.cd();c2.SetLogy(1);
  c2.cd();GGtframe->Draw();
  c2.SaveAs("SmearFactorfit.gif");

  return 1;
}
