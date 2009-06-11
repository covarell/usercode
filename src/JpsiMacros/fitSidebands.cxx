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
  //LEFT
  RooDataSet *GGleftside  = (RooDataSet*)data->reduce("JpsiMass < 2.9 && JpsiType == JpsiType::GG");
  GGleftside->setWeightVar(MCweight);
  RooRealVar GGcoefLeft("c_GG_Left","GG linear coefficient of left PDF",0.,-999.,999.);
  RooPolynomial GGleftFunct("GGleftFunc","GGleftFunc",JpsiMass,GGcoefLeft);
  GGleftFunct.fitTo(*GGleftside,Range("left"));

  RooPlot *GGleftmframe = JpsiMass.frame(Range("left"));
  GGleftmframe->SetTitle("Left sideband fit for glb-glb muons");
  GGleftside->plotOn(GGleftmframe,DataError(RooAbsData::SumW2));
  GGleftFunct.plotOn(GGleftmframe);
  GGleftFunct.paramOn(GGleftmframe);

  TCanvas c1;
  c1.cd();GGleftmframe->Draw();
  c1.SaveAs("GGbkmass_left.gif");

  //RIGHT
  RooDataSet *GGrightside = (RooDataSet*)data->reduce("JpsiMass > 3.3 && JpsiType == JpsiType::GG");
  GGrightside->setWeightVar(MCweight);
  RooRealVar GGcoefRight("c_GG_Right","GG linear coefficient of right PDF",0.,-999.,999.);
  RooPolynomial GGrightFunct("GGrightFunc","GGrightFunc",JpsiMass,GGcoefRight);
  GGrightFunct.fitTo(*GGrightside,Range("right"));

  RooPlot *GGrightmframe = JpsiMass.frame(Range("right"));
  GGrightmframe->SetTitle("Right sideband fit for glb-glb muons");
  GGrightside->plotOn(GGrightmframe,DataError(RooAbsData::SumW2));
  GGrightFunct.plotOn(GGrightmframe);
  GGrightFunct.paramOn(GGrightmframe);

  TCanvas c2;
  c2.cd();GGrightmframe->Draw();
  c2.SaveAs("GGbkmass_right.gif");

  //GT
  //LEFT
  RooDataSet *GTleftside  = (RooDataSet*)data->reduce("JpsiMass < 2.9 && JpsiType == JpsiType::GT");
  GTleftside->setWeightVar(MCweight);
  RooRealVar GTcoefLeft("c_GT_Left","GT linear coefficient of left PDF",0.,-999.,999.);
  RooPolynomial GTleftFunct("GTleftFunc","GTleftFunc",JpsiMass,GTcoefLeft);
  GTleftFunct.fitTo(*GTleftside,Range("left"));

  RooPlot *GTleftmframe = JpsiMass.frame(Range("left"));
  GTleftmframe->SetTitle("Left sideband fit for glb-trk muons");
  GTleftside->plotOn(GTleftmframe,DataError(RooAbsData::SumW2));
  GTleftFunct.plotOn(GTleftmframe);
  GTleftFunct.paramOn(GTleftmframe);

  TCanvas c3;
  c3.cd();GTleftmframe->Draw();
  c3.SaveAs("GTbkmass_left.gif");

  //RIGHT
  RooDataSet *GTrightside = (RooDataSet*)data->reduce("JpsiMass > 3.3 && JpsiType == JpsiType::GT");
  GTrightside->setWeightVar(MCweight);
  RooRealVar GTcoefRight("c_GT_Right","GT linear coefficient of right PDF",0.,-999.,999.);
  RooPolynomial GTrightFunct("GTrightFunc","GTrightFunc",JpsiMass,GTcoefRight);
  GTrightFunct.fitTo(*GTrightside,Range("right"));

  RooPlot *GTrightmframe = JpsiMass.frame(Range("right"));
  GTrightmframe->SetTitle("Right sideband fit for glb-trk muons");
  GTrightside->plotOn(GTrightmframe,DataError(RooAbsData::SumW2));
  GTrightFunct.plotOn(GTrightmframe);
  GTrightFunct.paramOn(GTrightmframe);

  TCanvas c4;
  c4.cd();GTrightmframe->Draw();
  c4.SaveAs("GTbkmass_right.gif");

  /// Global plots
  //GG
  RooPlot *GGmframe = JpsiMass.frame();
  GGmframe->SetTitle("#mu^{+}#mu^{-} mass for all glb-glb muons");
  data->plotOn(GGmframe,Cut("JpsiType == JpsiType::GG"),DataError(RooAbsData::SumW2));

  RooPlot *GGbackmframe = JpsiMass.frame(40);
  GGbackmframe->SetTitle("#mu^{+}#mu^{-} mass for background glb-glb muons");
  data->plotOn(GGbackmframe,Cut("JpsiType == JpsiType::GG && MCType == MCType::BK"),DataError(RooAbsData::SumW2));

  RooPlot *GGprmframe = JpsiMass.frame();
  GGprmframe->SetTitle("#mu^{+}#mu^{-} mass for prompt glb-glb muons");
  data->plotOn(GGprmframe,Cut("JpsiType == JpsiType::GG && MCType == MCType::PR"),DataError(RooAbsData::SumW2));

  RooPlot *GGnpmframe = JpsiMass.frame();
  GGprmframe->SetTitle("#mu^{+}#mu^{-} mass for non-prompt glb-glb muons");
  data->plotOn(GGnpmframe,Cut("JpsiType == JpsiType::GG && MCType == MCType::NP"),DataError(RooAbsData::SumW2));

  TCanvas c5;
  c5.Divide(2,2);
  c5.cd(1);GGmframe->Draw();
  c5.cd(2);GGbackmframe->Draw();
  c5.cd(3);GGprmframe->Draw();
  c5.cd(4);GGnpmframe->Draw();
  c5.SaveAs("GGbkmasses.gif");

  //GT
  RooPlot *GTmframe = JpsiMass.frame();
  GTmframe->SetTitle("#mu^{+}#mu^{-} mass for all glb-trk muons");
  data->plotOn(GTmframe,Cut("JpsiType == JpsiType::GT"),DataError(RooAbsData::SumW2));

  RooPlot *GTbackmframe = JpsiMass.frame();
  GTbackmframe->SetTitle("#mu^{+}#mu^{-} mass for background glb-trk muons");
  data->plotOn(GTbackmframe,Cut("JpsiType == JpsiType::GT && MCType == MCType::BK"),DataError(RooAbsData::SumW2));

  RooPlot *GTprmframe = JpsiMass.frame();
  GTprmframe->SetTitle("#mu^{+}#mu^{-} mass for prompt glb-trk muons");
  data->plotOn(GTprmframe,Cut("JpsiType == JpsiType::GT && MCType == MCType::PR"),DataError(RooAbsData::SumW2));

  RooPlot *GTnpmframe = JpsiMass.frame();
  GTprmframe->SetTitle("#mu^{+}#mu^{-} mass for non-prompt glb-trk muons");
  data->plotOn(GTnpmframe,Cut("JpsiType == JpsiType::GT && MCType == MCType::NP"),DataError(RooAbsData::SumW2));

  TCanvas c6;
  c6.Divide(2,2);
  c6.cd(1);GTmframe->Draw();
  c6.cd(2);GTbackmframe->Draw();
  c6.cd(3);GTprmframe->Draw();
  c6.cd(4);GTnpmframe->Draw();
  c6.SaveAs("GTbkmasses.gif");

  return 1;
}
