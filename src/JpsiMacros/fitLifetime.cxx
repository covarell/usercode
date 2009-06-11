// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooCategory.h>
#include <RooPolynomial.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooConstVar.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooHist.h>

using namespace::RooFit;

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

  RooRealVar Jpsict("Jpsict","#mu^{+}#mu^{-} c#tau",JpsiCtMin,JpsiCtMax,"mm");

  RooRealVar MCweight("MCweight","Monte Carlo Weight",0.,5.);

  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  //GT
  //Left sideband
  RooDataSet *GTleftdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT && JpsiMass < 2.9");
  RooDataSet *GTrightdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GT && JpsiMass > 3.2");
  GTleftdata->setWeightVar(MCweight);
  GTrightdata->setWeightVar(MCweight);

  RooDataSet *GGleftdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG && JpsiMass < 2.9");
  RooDataSet *GGrightdata = (RooDataSet*)data->reduce("JpsiType == JpsiType::GG && JpsiMass > 3.2");
  GGleftdata->setWeightVar(MCweight);
  GGrightdata->setWeightVar(MCweight);

  RooRealVar meanSig("meanSig","Mean of the resolution gaussian",0.,-1.,1.);
  RooRealVar sigmaSig("sigmaSig","#sigma of the resolution gaussian",0.1,0.,5.);
  RooRealVar scaleK("scaleK","Scale factor of the resolution gaussian",1.,0.,10.);
  RooConstVar one("one","one",1.);

  scaleK.setConstant(kTRUE);

  RooRealVar lambdap("lambdap","tau of the positive background tail",0.1,0.,5.);
  RooRealVar lambdam("lambdam","tau of the negative background tail",0.1,0.,5.);
  RooRealVar lambdasym("lambdasym","tau of the symmetric background tail",0.1,0.,5.);

  RooGaussModel resG("resG","Gaussian resolution function",Jpsict,meanSig,sigmaSig,one,scaleK);

  RooFormulaVar denomGauss("denomGauss","@0*@1",RooArgList(sigmaSig,scaleK));

  RooGaussian bkg1("bkg1","Zero lifetime background",Jpsict,meanSig,denomGauss);
  RooDecay bkg2("bkg2","One sided positive background",Jpsict,lambdap,resG,RooDecay::SingleSided);
  RooDecay bkg3("bkg3","One sided negative background",Jpsict,lambdam,resG,RooDecay::Flipped);
  RooDecay bkg4("bkg4","Symmetric background",Jpsict,lambdasym,resG,RooDecay::DoubleSided);

  //now, we could just compose the background together
  //but since we are not idiots, we compose them partially for stability
  RooRealVar fpm("fpm","Fraction of pos/neg tails",0.5,0.,1.);
  RooAddPdf bkgPart1("bkgPart1","Sum of pos/neg backgrounds",bkg2,bkg3,fpm);

  RooRealVar fLiving("fLiving","Fraction of sym/asym living backgrounds",0.5,0.,1.);
  RooAddPdf bkgPart2("bkgPart2","Sum of living backgrounds",bkgPart1,bkg4,fLiving);

  RooRealVar fTot("fTot","Fraction of delta living background",0.5,0.,1.);
  RooAddPdf bkgTOT("bkgTOT","Sum of all backgrounds",bkg1,bkgPart2,fTot);

  //now fit the left m_mumu sideband in ctau

  bkgTOT.fitTo(*GTleftdata,Minos(0));

  RooPlot *GTtleftframe = Jpsict.frame();
  GTtleftframe->SetTitle("c#tau fit for glb-trk muons in left mass sideband");
  GTleftdata->plotOn(GTtleftframe,DataError(RooAbsData::SumW2));
  bkgTOT.plotOn(GTtleftframe);

  RooPlot* pullGTleftT = Jpsict.frame() ;
  pullGTleftT->addPlotable(GTtleftframe->pullHist()) ;
  pullGTleftT->SetMaximum(5.);
  pullGTleftT->SetMinimum(-5.);
  pullGTleftT->SetTitle("Residuals");

  TCanvas *c1 = new TCanvas("c1","c1",400,600);

  TPad *c1_1 = new TPad("c1_1","c1_1",0.,0.3,1.,1.);
  TPad *c1_2 = new TPad("c1_2","c1_2",0.,0.,1.,0.3);

  c1_1->Draw();
  c1_2->Draw();

  c1_1->cd();c1_1->SetLogy(1);
  GTtleftframe->Draw();

  c1_2->cd();pullGTleftT->Draw();

  c1->SaveAs("GTbkgctFit_left.gif");

  bkgTOT.fitTo(*GTrightdata,Minos(0));

  RooPlot *GTtrightframe = Jpsict.frame();
  GTtrightframe->SetTitle("c#tau fit for glb-trk muons in right mass sideband");
  GTrightdata->plotOn(GTtrightframe,DataError(RooAbsData::SumW2));
  bkgTOT.plotOn(GTtrightframe);

  RooPlot* pullGTrightT = Jpsict.frame() ;
  pullGTrightT->addPlotable(GTtrightframe->pullHist()) ;
  pullGTrightT->SetMaximum(5.);
  pullGTrightT->SetMinimum(-5.);
  pullGTrightT->SetTitle("Residuals");

  TCanvas *c2 = new TCanvas("c1","c1",400,600);

  TPad *c2_1 = new TPad("c2_1","c2_1",0.,0.3,1.,1.);
  TPad *c2_2 = new TPad("c2_2","c2_2",0.,0.,1.,0.3);

  c2_1->Draw();
  c2_2->Draw();

  c2_1->cd();c2_1->SetLogy(1);
  GTtrightframe->Draw();

  c2_2->cd();pullGTrightT->Draw();

  c2->SaveAs("GTbkgctFit_right.gif");

  bkgTOT.fitTo(*GGleftdata,Minos(0));

  RooPlot *GGtleftframe = Jpsict.frame();
  GGtleftframe->SetTitle("c#tau fit for glb-glb muons in left mass sideband");
  GGleftdata->plotOn(GGtleftframe,DataError(RooAbsData::SumW2));
  bkgTOT.plotOn(GGtleftframe);

  TCanvas c3;
  c3.cd();c3.SetLogy(1);
  GGtleftframe->Draw();
  c3.SaveAs("GGbkgctFit_left.gif");

  bkgTOT.fitTo(*GGrightdata,Minos(0));

  RooPlot *GGtrightframe = Jpsict.frame();
  GGtrightframe->SetTitle("c#tau fit for glb-glb muons in right mass sideband");
  GGrightdata->plotOn(GGtrightframe,DataError(RooAbsData::SumW2));
  bkgTOT.plotOn(GGtrightframe);

  TCanvas c4;
  c4.cd();c4.SetLogy(1);
  GGtrightframe->Draw();
  c4.SaveAs("GGbkgctFit_right.gif");


  return 1;
}
