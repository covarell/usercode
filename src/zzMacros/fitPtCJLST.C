//Fit the PT curves for sig and BKG
//
//By: Roberto Covarelli (U of Rochester)

#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TPad.h"
#include "TLatex.h"

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooBinning.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooRandom.h"

// #include "PDFs/RooTsallis.h"
// #include "PDFs/RooTsallis2.h"
#include "PDFs/RooTsallis3.h" 

using namespace RooFit;

RooDataHist* withSmartBinning(TH1F* source, RooRealVar* var, float min, float max, float minWeight, float sizeSig, int type) {
  
  // DEFINE BINNING
  // type:
  // 0 = background
  // 1 = HRes gg signal
  // 2 = POWHEG gg signal
  // 2 = VBF signal
  float lowlimit;
  std::vector<float> sizes;
  std::vector<float> highlimits;  
  if (type == 0) {
    lowlimit = min;
    sizes.push_back(1.);
    highlimits.push_back(max);
  } else if (type == 1 || type == 2) {
    lowlimit = min;
    sizes.push_back(sizeSig);        
    // sizes.push_back(8.);        
    // sizes.push_back(4.);
    // sizes.push_back(8.);
    // highlimits.push_back(150.); 
    // highlimits.push_back(200.);
    // highlimits.push_back(320.);
    highlimits.push_back(max);
  } else {
    lowlimit = min;
    sizes.push_back(sizeSig);
    highlimits.push_back(max);
  }

  // ROOREALVAR BINNING
  var->setMin(lowlimit);
  var->setMax(highlimits.back());
  RooBinning rb(lowlimit,highlimits.back());
  float leftlimit = lowlimit;
  for (unsigned int ii = 0; ii < highlimits.size(); ii++) {
    float rightlimit = highlimits.at(ii);
    int howManyBins = (rightlimit-leftlimit)/sizes.at(ii);
    rb.addUniform(howManyBins,leftlimit,rightlimit);
    leftlimit = rightlimit;
  }
  var->setBinning(rb);

  // FILL ROODATAHIST
  RooDataHist* result = new RooDataHist("result","A dataset",RooArgList(*var));
  float nWeight = source->GetEntries();
  if (type == 2) nWeight = 0.00215;
  float binWidth = source->GetBinWidth(1);

  float thisWeight = 0.;
  float thisWeightsum = 0.;
  int j = 0;
  int whichInterval = 0;
  for (Int_t i=1; i<=source->GetNbinsX(); i++) {

    float thispt = source->GetXaxis()->GetBinCenter(i);
    if (thispt > lowlimit && thispt < highlimits.back()) {
      bool lastbin = false;
      if (source->GetXaxis()->GetBinCenter(i+1) >= highlimits.back()) lastbin = true;
      if (thispt > highlimits.at(whichInterval)) whichInterval++;
      var->setVal(thispt);
      thisWeight += source->GetBinContent(i)*nWeight;
      thisWeightsum += source->GetBinError(i)*nWeight*source->GetBinError(i)*nWeight;
      j++;
      cout << "Entries in TH1 bin " << i << " = " << source->GetBinContent(i)*nWeight << " +/- " << source->GetBinError(i)*nWeight << endl;   
      if (binWidth*j >= sizes.at(whichInterval) || lastbin) {
        // cout << "Entro" << endl; 
	if (thisWeight <= 0) thisWeight = minWeight; 
        if (fabs(thisWeight/sizes.at(whichInterval) - 0.0829038) < 0.00001) thisWeight /= 10.;
	if (sqrt(thisWeightsum) > thisWeight) thisWeightsum = thisWeight*thisWeight*0.9;
	cout << " " << var->getVal() << " " << thisWeight << " " << thisWeightsum << " " << sizes.at(whichInterval) << endl; 
	result->add(RooArgSet(*var),
		    thisWeight/sizes.at(whichInterval),
		    thisWeightsum/(sizes.at(whichInterval)*sizes.at(whichInterval)));
	thisWeight = 0.;
	thisWeightsum = 0.;
	j = 0;
      }
    }
  } 

  return result;
}

void fitPtCJLST(int LHCsqrts = 7, int whichtype = 1, 
		bool correctErrors = false, string systString = "Default")

// whichtype
// 0 - gg Signal
// 1 - VBF Signal
// 2 - ZZ
// 3 = ZX

// So far only for 125 GeV...

{

  string nameSample[4] = {"gg","vbf","zz","zx"};
  float maxType[4] = {300.,400.,250.,200.};   // GeV
  float rebinType[4] = {1,2,1,1};
  
  char fileToOpen[200];
  sprintf(fileToOpen,"PT_%s_SEL_%dTeV.root",nameSample[whichtype].c_str(),LHCsqrts);
  if (whichtype == 3) sprintf(fileToOpen,"PT_%s_SEL_allTeV.root",nameSample[whichtype].c_str());

  RooRealVar* pt = new RooRealVar("pt","p_{T}^{H}",0.,maxType[whichtype],"GeV/c");
 
  TFile input(fileToOpen);
  sprintf(fileToOpen,"ptH_%s",systString.c_str());
  TH1F* ptH = (TH1F*)input.Get(fileToOpen);
  
  if (rebinType[whichtype] > 1) ptH->Rebin(rebinType[whichtype]);
  if (maxType[whichtype] < ptH->GetBinLowEdge(ptH->GetNbinsX() + 1) - ptH->GetBinWidth(1)) {
    int theBin = ptH->FindBin(maxType[whichtype]);
    ptH->GetXaxis()->SetRange(1,theBin-1);
  }

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // cout << endl << "Signal " << endl;   
  pt->setBins(ptH->GetNbinsX());

  RooDataHist* rdh = new RooDataHist("rdh","Some dataset",RooArgList(*pt),Import(*ptH,kFALSE));
 
  // fit definitions
  RooWorkspace *ws = new RooWorkspace("ws");

  RooRealVar m("m","emme", 110.,10., 1200.);
  RooRealVar n("n","enne", 0.93, 0.5, 15.);
  RooRealVar n2("n2","enne2", 0.75, 0.5, 15.);
  RooRealVar bb("bb","bibi",0.02, 0.0005, 0.1);
  RooRealVar T("T","tti",0.2,0.000005,1.);
  RooRealVar bb2("bb2","bibi2",0.02, 0.0005, 0.1);
  RooRealVar fexp("fexp","f_exp",0.02, 0.0, 1.0);
  if (whichtype == 2) {
    if (LHCsqrts == 8) {
      m.setVal(79.885);   // m.setConstant(kTRUE);    
      bb.setVal(0.020); // bb.setConstant(kTRUE);
      n2.setVal(1.0678);  // n2.setConstant(kTRUE);
      n.setVal(1.010);    // n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.088077);   T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    } else {
      m.setVal(83.54);   // m.setConstant(kTRUE);    
      bb.setVal(0.020); // bb.setConstant(kTRUE);
      n2.setVal(1.0803);  // n2.setConstant(kTRUE);
      n.setVal(1.0207);   //  n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.22341);    T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    }
  }
  else if (whichtype == 1) {
    m.setVal(235.3);   // m.setConstant(kTRUE);
    n.setVal(1.1705);   // n.setConstant(kTRUE);
    n2.setVal(4.086);   n2.setConstant(kTRUE);
    bb.setVal(0.0294); // bb.setConstant(kTRUE);
    T.setVal(0.0000703);    T.setConstant(kTRUE);
    bb2.setVal(0.0203);   bb2.setConstant(kTRUE);
    fexp.setVal(0.100);   fexp.setConstant(kTRUE);
  }
  else if (whichtype == 3) {
    m.setVal(653.8);   // m.setConstant(kTRUE);
    n.setVal(0.8367);   n.setConstant(kTRUE);
    n2.setVal(2.985);   n2.setConstant(kTRUE);
    bb.setVal(0.0562); // bb.setConstant(kTRUE);
    T.setVal(0.0000938);   // T.setConstant(kTRUE);
    bb2.setVal(0.01037);   bb2.setConstant(kTRUE);
    fexp.setVal(0.100);   fexp.setConstant(kTRUE);
  } else {
    if (LHCsqrts == 8) {
      m.setVal(591.06);   // m.setConstant(kTRUE);
      n.setVal(0.6048);   // if (systString == "Default") n.setConstant(kTRUE);
      n2.setVal(0.9086);   n2.setConstant(kTRUE);
      bb.setVal(0.0280);   // bb.setConstant(kTRUE);
      T.setVal(0.0866);   // T.setConstant(kTRUE);
      bb2.setVal(0.00657);   bb2.setConstant(kTRUE);
      fexp.setVal(0.0849);   // if (systString == "Default") fexp.setConstant(kTRUE);
    } else {
      if (systString == "Resummation") m.setMax(800.);
      m.setVal(456.99);   // m.setConstant(kTRUE);
      n.setVal(1.068);   // n.setConstant(kTRUE);
      n2.setVal(0.9649);   n2.setConstant(kTRUE);
      bb.setVal(0.0149);   // bb.setConstant(kTRUE);
      T.setVal(0.2330);   // T.setConstant(kTRUE);
      bb2.setVal(0.00253);    bb2.setConstant(kTRUE);
      fexp.setVal(0.0202);   fexp.setConstant(kTRUE);
    }
  }
  
  RooTsallis3* rt3 = new RooTsallis3("rt3","rt3",*pt,m,n,n2,bb,bb2,T,fexp);
  ws->import(*rt3);

  // fit
  RooFitResult* fit = ws->pdf("rt3")->fitTo(*rdh,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  RooPlot *frame = pt->frame();

  char reducestr[300];
  sprintf(reducestr,"pt > %f && pt < %f",pt->getMin(),pt->getMax());
  
  rdh->plotOn(frame,DataError(RooAbsData::SumW2),Cut(reducestr));
  ws->pdf("rt3")->plotOn(frame,LineColor(kBlue),Normalization(rdh->sumEntries(),RooAbsReal::NumEvent));
  RooHist *hpull = frame->pullHist();
  float chi2 = 0.;

  double *ypulls = hpull->GetY();
  unsigned int nBins = rdh->numEntries();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    cout << "Pull of bin " << i << " = " << ypulls[i] << endl;
    if (fabs(ypulls[i]) < 5.0) chi2 += ypulls[i]*ypulls[i]; 
    cout << "Partial chi2 = " << chi2 << endl;
    if (fabs(ypulls[i]) > 0.0001 && fabs(ypulls[i]) < 5.0) nFullBins++;
  }
  for (unsigned int i = 0; i < nBins; i++) {
    if (fabs(ypulls[i]) < 0.0001) ypulls[i] = 999.; 
    hpull->SetPointError(i,0.,0.,0.,0.);
  } 
  int nFitPar = fit->floatParsFinal().getSize() - 1;

  if (correctErrors) {
    // Tsallis errors not reliable, use toy MC

    float mVal = m.getVal();
    float nVal = n.getVal();
    float n2Val = n2.getVal();
    float bbVal = bb.getVal();
    float bb2Val = bb2.getVal();
    float fexpVal = fexp.getVal();
    float TVal = T.getVal();

    TH1F* mHist = new TH1F("mHist","m",21,-0.5*mVal,0.5*mVal);
    TH1F* nHist = new TH1F("nHist","n",21,-0.2*nVal,0.2*nVal);
    TH1F* n2Hist = new TH1F("n2Hist","n2",21,-0.2*n2Val,0.2*n2Val);
    TH1F* bbHist = new TH1F("bbHist","bb",21,-0.2*bbVal,0.2*bbVal);
    TH1F* bb2Hist = new TH1F("bb2Hist","bb2",21,-0.2*bb2Val,0.2*bb2Val);
    TH1F* fexpHist = new TH1F("fexpHist","fexp",21,-0.2*fexpVal-0.000001,0.2*fexpVal+0.000001);
    TH1F* THist = new TH1F("THist","T",21,-0.5*TVal,0.5*TVal);

    for (unsigned int iToy = 0; iToy < 200; iToy++) {

      cout << endl << "####" << endl;
      cout << "Generating toy experiment n. " << iToy+1 << endl;

      m.setVal(mVal);
      n.setVal(nVal);
      n2.setVal(n2Val);
      bb.setVal(bbVal);
      bb2.setVal(bb2Val);
      fexp.setVal(fexpVal);
      T.setVal(TVal);

      TDatime *now = new TDatime();
      Int_t seed = now->GetDate() + now->GetTime();
      cout << "RooFit Generation Seed = " << seed+iToy << endl;
      RooRandom::randomGenerator()->SetSeed(seed+iToy);
      cout << "####" << endl << endl;

      RooDataSet *dataToy = rt3->generate(RooArgSet(*pt),ptH->GetEntries());
      RooDataHist *dataToyH = new RooDataHist("dataToyH","toy",RooArgSet(*pt),*dataToy);
      
      rt3->fitTo(*dataToyH,Minos(0),SumW2Error(kTRUE),NumCPU(1));  
  
      if (fit->floatParsFinal().find("m")) mHist->Fill(m.getVal()-mVal);
      if (fit->floatParsFinal().find("n")) nHist->Fill(n.getVal()-nVal);
      if (fit->floatParsFinal().find("n2")) n2Hist->Fill(n2.getVal()-n2Val);
      if (fit->floatParsFinal().find("bb")) bbHist->Fill(bb.getVal()-bbVal);
      if (fit->floatParsFinal().find("bb2")) bbHist->Fill(bb2.getVal()-bb2Val);
      if (fit->floatParsFinal().find("fexp")) fexpHist->Fill(fexp.getVal()-fexpVal);
      if (fit->floatParsFinal().find("T")) THist->Fill(T.getVal()-TVal);
    }

    TCanvas cant("cant","Test canvas",5.,5.,900.,500.);
    cant.Divide(4,2);
    cant.cd(1);   mHist->Draw();
    cant.cd(2);   nHist->Draw();
    cant.cd(3);   n2Hist->Draw();
    cant.cd(4);   bbHist->Draw();
    cant.cd(5);   bb2Hist->Draw();
    cant.cd(6);   fexpHist->Draw();
    cant.cd(7);   THist->Draw();
    cant.SaveAs("figs/testToys.pdf");

    if (fit->floatParsFinal().find("m")) m.setError(mHist->GetRMS());
    if (fit->floatParsFinal().find("n")) n.setError(nHist->GetRMS());
    if (fit->floatParsFinal().find("n2")) n2.setError(n2Hist->GetRMS());
    if (fit->floatParsFinal().find("bb")) bb.setError(bbHist->GetRMS());
    if (fit->floatParsFinal().find("bb2")) bb2.setError(bb2Hist->GetRMS());
    if (fit->floatParsFinal().find("fexp")) fexp.setError(fexpHist->GetRMS());
    if (fit->floatParsFinal().find("T")) T.setError(THist->GetRMS());
  }

  char fileToSave[200];
  sprintf(fileToSave,"text/paramsCJLST_%s_%dTeV_%s.txt",nameSample[whichtype].c_str(),LHCsqrts,systString.c_str());
  ofstream os1(fileToSave);
  (RooArgSet(m,n,n2,bb,bb2,fexp,T)).writeToStream(os1,false);
  os1.close();

  TCanvas can("can","The canvas",5.,5.,500.,900.); 
  can.Divide(1,3);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);

  can.cd(1);
  gPad->SetBottomMargin(0.0);
  frame->Draw();
  // gPad->SetLogy(); 
  // Htest->Draw();
  sprintf(fileToSave,"%s - %d TeV",nameSample[whichtype].c_str(),LHCsqrts);
  t->DrawLatex(0.6,0.8,fileToSave); 

  can.cd(2);
  gPad->SetLogy(); 
  gPad->SetTopMargin(0.0);
  frame->Draw();
 
  RooPlot* pull = pt->frame(Title("Pull Distribution")) ;
  pull->GetYaxis()->SetTitle("Pull");
  /* pull->SetLabelSize(0.08,"XYZ");
  pull->SetTitleSize(0.08,"XYZ");
  pull->SetTitleOffset(0.6,"Y");
  pull->SetTitleOffset(1.0,"X"); */
  pull->addPlotable(hpull,"P") ; 
  pull->SetMinimum(-6.); 
  pull->SetMaximum(6.); 

  can.cd(3);
  gPad->SetGridy();
  pull->Draw();
  sprintf(fileToSave,"#chi^{2}/n_{DoF} = %4.1f/%d",chi2,nFullBins - nFitPar);
  if (chi2 < 1000.) t->DrawLatex(0.80,0.86,fileToSave);

  sprintf(fileToSave,"figs/fitCJLST_%s_%dTeV_%s.pdf",nameSample[whichtype].c_str(),LHCsqrts,systString.c_str());
  can.SaveAs(fileToSave);

}



