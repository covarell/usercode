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
#include "RooExponential.h"
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

using namespace RooFit;

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  bool useSidebands;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){
      case 'f':
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      case 's':
        useSidebands = true;
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
  const float JpsiCtMax = 3.0;

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

  //add the parameters to save/retrieve
  RooArgSet paramlist;

  //CONSIDER THE CASE
  RooDataSet *GTdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT");
  GTdata->setWeightVar(MCweight);
  RooDataSet *GTdataPr = (RooDataSet*)GTdata->reduce("MCType == MCType::PR");
  RooDataSet *GTdataNp = (RooDataSet*)GTdata->reduce("MCType == MCType::NP");
  GTdataPr->setWeightVar(MCweight);
  GTdataNp->setWeightVar(MCweight); 

  cout << "Number of events to fit  = " << GTdata->numEntries(kTRUE) << endl; 

  //JPSI MASS PARAMETRIZATION

  //background
  RooRealVar coefExp("coefExp","exponential coefficient of bkg PDF",-5.,-9.,0.);
  RooExponential GTexpFunct("GTexpFunc","GTexpFunc",JpsiMass,coefExp); 

  paramlist.add(coefExp);

  //signal is the same for prompt and non-prompt
  RooRealVar meanMassSig1("meanMassSig1","Mean of the signal gaussian 1",3.1,2.8,3.2);
  RooRealVar sigmaMassSig1("sigmaMassSig1","#sigma of the signal gaussian 1",0.02,0.,0.9);

  paramlist.add(meanMassSig1);
  paramlist.add(sigmaMassSig1);

  RooRealVar meanMassSig2("meanMassSig2","Mean of the signal gaussian 2",3.0,2.8,3.2);
  RooRealVar sigmaMassSig2("sigmaMassSig2","#sigma of the signal gaussian 2",0.04,0.,0.9);

  paramlist.add(meanMassSig2);
  paramlist.add(sigmaMassSig2);

  RooGaussian signalMassG1("signalMassG1","Signal PDF 1",JpsiMass,meanMassSig1,sigmaMassSig1);
  RooGaussian signalMassG2("signalMassG2","Signal PDF 2",JpsiMass,meanMassSig1,sigmaMassSig2);

  RooRealVar coeffGauss("coeffGauss","Relative norm of the two signal gaussians",0.36,0.,1.);

  paramlist.add(coeffGauss);

  RooAddPdf sigMassPDF("sigMadssPDF","Total signal pdf",signalMassG1,signalMassG2,coeffGauss);

  //JPSI CTAU PARAMETRIZATION

  //resolution function
  RooRealVar meanResSigW("meanResSigW","Mean of the resolution wide gaussian",0.,-1.,1.);
  RooRealVar sigmaResSigW("sigmaResSigW","#sigma of the resolution wide gaussian",0.6,0.,5.);
  RooRealVar scaleK("scaleK","Scale factor of the resolution gaussian",1.,0.,10.);
  RooConstVar one("one","one",1.);

  paramlist.add(meanResSigW);
  paramlist.add(sigmaResSigW);
  paramlist.add(scaleK);

  scaleK.setConstant(kTRUE);

  RooRealVar meanResSigN("meanResSigN","Mean of the resolution narrow gaussian",0.,-1.,1.);
  RooRealVar sigmaResSigN("sigmaResSigN","#sigma of the resolution narrow gaussian",0.1,0.,5.);

  paramlist.add(meanResSigN);
  paramlist.add(sigmaResSigN);

  RooGaussModel resGW("resGW","Wide Gaussian resolution function",Jpsict,meanResSigW,sigmaResSigW,one,scaleK);
  RooGaussModel resGN("resGN","Narrow Gaussian resolution function",Jpsict,meanResSigN,sigmaResSigN,one,scaleK);

  RooRealVar fracRes("fracRes","Fraction of narrow/wider gaussians",0.5,0.,1.);

  paramlist.add(fracRes);

  RooAddModel resol("resol","resol",RooArgList(resGW,resGN),RooArgList(fracRes));

  //background
  RooRealVar lambdap("lambdap","tau of the positive background tail",0.40,0.,5.);
  RooRealVar lambdam("lambdam","tau of the negative background tail",0.12,0.,5.);
  RooRealVar lambdasym("lambdasym","tau of the symmetric background tail",1.37,0.,5.);

  paramlist.add(lambdap);
  paramlist.add(lambdam);
  paramlist.add(lambdasym);

  //bkg1 is the resolution function

  RooDecay bkg2("bkg2","One sided positive background",Jpsict,lambdap,resol,RooDecay::SingleSided);
  RooDecay bkg3("bkg3","One sided negative background",Jpsict,lambdam,resol,RooDecay::Flipped);
  RooDecay bkg4("bkg4","Symmetric background",Jpsict,lambdasym,resol,RooDecay::DoubleSided);

  //now, we could just compose the background together
  //but since we are not idiots, we compose them partially for stability
  RooRealVar fpm("fpm","Fraction of pos/neg tails",0.68,0.,1.);

  paramlist.add(fpm);

  RooAddPdf bkgPart1("bkgPart1","Sum of pos/neg backgrounds",bkg2,bkg3,fpm);

  RooRealVar fLiving("fLiving","Fraction of sym/asym living backgrounds",0.94,0.,1.);

  paramlist.add(fLiving);

  RooAddPdf bkgPart2("bkgPart2","Sum of living backgrounds",bkgPart1,bkg4,fLiving);

  RooRealVar fbkgTot("fbkgTot","Fraction of delta living background",0.08,0.,1.);

  paramlist.add(fbkgTot);

  RooAddPdf bkgctauTOT("bkgctauTOT","Sum of all backgrounds",resol,bkgPart2,fbkgTot);

  //signal prompt, same as zero lifetime background

  //signal non-prompt
  RooRealVar taueff("taueff","Effective tau of the B meson",0.348,0.,1.);
  //RooRealVar sigmaMC("sigmaMC","MC #sigma",0.3,0.,5.);

  paramlist.add(taueff);
  //paramlist.add(sigmaMC);

  //RooGExpModel physsigNP("physsigNP","Gauss + exp model",Jpsict,sigmaMC,taueff);

  //RooUniformBinning FFTbin(Jpsict.getMin(),Jpsict.getMax(),1000,"cache");
  //Jpsict.setBinning(FFTbin,"cache");
  //RooFFTConvPdf sigNP("sigNP","Non-prompt signal",Jpsict,physsigNP,resol);

  RooDecay sigNP("sigNP","Non-prompt signal",Jpsict,taueff,resol,RooDecay::SingleSided);

  //putting all together
  RooProdPdf totsigPR("totsigPR","Total prompt signal",RooArgList(sigMassPDF,resol));
  RooProdPdf totsigNP("totsigNP","Total non-prompt signal",RooArgList(sigMassPDF,sigNP));
  RooProdPdf totBKG("totBKG","Total background",RooArgList(GTexpFunct,bkgctauTOT));

  RooRealVar NSigPR("NSigPR","Number of prompt signal events",4000.,10.,1000000.);
  RooRealVar NSigNP("NSigNP","Number of non-prompt signal events",900.,10.,1000000.);
  RooRealVar NBkg("NBkg","Number of background events",1400.,10.,1000000.);

  paramlist.add(NSigPR);
  paramlist.add(NSigNP);
  paramlist.add(NBkg);

  RooAddPdf totPDF("totPDF","Total PDF",RooArgList(totsigPR,totsigNP,totBKG),RooArgList(NSigPR,NSigNP,NBkg));

  paramlist.readFromFile("fit2dpars_GT.txt");

  RooFitResult* fitRes = totPDF.fitTo(*GTdata,Extended(1),Save(1),Minos(0));

  Double_t Bfrac = NSigNP.getVal()/(NSigNP.getVal() + NSigPR.getVal());
  cout << "B frac = " << Bfrac << " +/- " << Bfrac/NSigNP.getVal() << endl;

  RooArgSet results(fitRes->floatParsFinal());
  RooArgSet conresults(fitRes->constPars());
  results.add(conresults);
  results.writeToFile("fit2d_results_GT.txt");

  RooPlot *GTmframe = JpsiMass.frame();
  GTmframe->SetTitle("2D fit for glb-trk muons (mass projection)");
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GTmframe);
  totPDF.plotOn(GTmframe,DrawOption("F"),FillColor(kGreen));
  totPDF.plotOn(GTmframe,Components(RooArgSet(totsigNP,totBKG)),DrawOption("F"),FillColor(kBlue));
  totPDF.plotOn(GTmframe,Components(RooArgSet(totBKG)),DrawOption("F"),FillColor(kRed));
  GTdata->plotOn(GTmframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GTmframe);

  TCanvas c1;
  c1.cd();GTmframe->Draw();
  c1.SaveAs("2DGTmassfit.gif");

  RooPlot *GTtframe = Jpsict.frame();
  GTtframe->SetTitle("2D fit for glb-trk muons (c   #tau projection)");
  GTdata->plotOn(GTtframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GTtframe);
  totPDF.plotOn(GTtframe,DrawOption("F"),FillColor(kGreen));
  totPDF.plotOn(GTtframe,Components(RooArgSet(totsigNP,totBKG)),DrawOption("F"),FillColor(kBlue));
  totPDF.plotOn(GTtframe,Components(RooArgSet(totBKG)),DrawOption("F"),FillColor(kRed));
  GTdata->plotOn(GTtframe,DataError(RooAbsData::SumW2));
  totPDF.plotOn(GTtframe);

  TCanvas c2;
  c2.cd(); c2.SetLogy(1);
  c2.cd();GTtframe->Draw();
  c2.SaveAs("2DGTtimefit.gif");

  cout << endl << "GT J/psi yields:" << endl;
  cout << "PROMPT :     True MC : " << GTdataPr->numEntries(true) << " Fit : " << NSigPR.getVal() << " +/- " << NSigPR.getError() << endl;
  cout << "NON-PROMPT : True MC : " << GTdataNp->numEntries(true) << " Fit : " << NSigNP.getVal() << " +/- " << NSigNP.getError() << endl;

  return 1;
}
