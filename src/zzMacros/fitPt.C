//Fit the PT curves for sig and BKG
//
//By: Roberto Covarelli (U of Rochester)

#include <math.h>
#include <string>
#include <iostream>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooFitResult.h"

// #include "PDFs/RooTsallis.h"
#include "PDFs/RooTsallis2.h"
#include "PDFs/RooTsallis3.h"

using namespace RooFit;

void fitPt(float mZZcenter = 126., float mZZspread = 5.)
{

  // range limits
  float varMin = 1.;
  float varMax = 120.;
  float binWidthBkg = 1.;
  float binWidthSig = 2.;
  float scaleErrorSig = 1.;  /// ???

  TFile* file = new TFile("PT_Y_Temp.root");
  TH1F* massH = (TH1F*)((TH2F*)file->Get("Pt_bkg"))->ProjectionX();
  float mZZmin = mZZcenter - mZZspread;
  float mZZmax = mZZcenter + mZZspread;
  int binMin = massH->FindBin(mZZmin);
  int binMax = massH->FindBin(mZZmax);

  TH1F* bkgH = (TH1F*)((TH2F*)file->Get("Pt_bkg"))->ProjectionY("bkgH",binMin,binMax);
  // TH1F* sigH = (TH1F*)((TH2F*)file->Get("Pt_sig"))->ProjectionY("sigH",binMin,binMax);
  
  char fileToOpen[200];
  sprintf(fileToOpen,"PT_%d_Temp.root",int(mZZcenter));
  TFile* file2 = new TFile(fileToOpen);
  TH1F* sigH = (TH1F*)file2->Get("Pt_sig");

  float nSigWeight = sigH->GetEntries(); 
 
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // Build datasets
  RooRealVar pt("pt","p_{T}^{H}",varMin, varMax,"GeV/c");
  int nbins = int((varMax-varMin)/binWidthBkg);
  pt.setBins(nbins);

  RooDataHist* bkg = new RooDataHist("bkg","Background dataset",RooArgList(pt),bkgH);

  int nbins2 = int((varMax-varMin)/binWidthSig);
  pt.setBins(nbins2);

  RooDataHist* sig = new RooDataHist("sig","Signal dataset",RooArgList(pt));
  
  for (Int_t i=1; i<=sigH->GetNbinsX(); i++) {

    float thispt = sigH->GetXaxis()->GetBinCenter(i);
    if (thispt > varMin && thispt < varMax) {
      pt.setVal(thispt);
      cout << "Entries in bin " << i << " = " << sigH->GetBinContent(i)*nSigWeight << endl;
      float weightsum = sigH->GetBinError(i)*nSigWeight*scaleErrorSig;
      // avoid negative entries
      if (weightsum > sigH->GetBinContent(i)*nSigWeight) weightsum = sigH->GetBinContent(i)*nSigWeight - 1;
      sig->add(RooArgSet(pt),sigH->GetBinContent(i)*nSigWeight,weightsum*weightsum);
    }
  }

  // fit definitions
  RooWorkspace *ws = new RooWorkspace("ws");

  RooRealVar m("m","emme", 110.,10., 600.,"GeV/c");
  RooRealVar n("n","enne", 0.93, 0.5, 15.);
  RooRealVar n2("n2","enne2", 0.75, 0.5, 15.);
  RooRealVar bb("bb","bibi",0.02, 0.0005, 0.1);
  RooRealVar T("T","tti",0.2,0.0005,1.);
  m.setVal(54.47);   m.setConstant(kTRUE);
  n.setVal(1.255);   // n.setConstant(kTRUE);
  n2.setVal(1.73);   n2.setConstant(kTRUE);
  bb.setVal(0.020); // bb.setConstant(kTRUE);
  T.setVal(0.20);   // T.setConstant(kTRUE);

  RooRealVar ms("ms","emme signal", 110.,10., 6000.,"GeV/c");
  RooRealVar ns("ns","enne signal", 0.93, 0.5, 15.);
  RooRealVar n2s("n2s","enne2 signal", 0.75, 0.5, 15.); 
  RooRealVar bbs("bbs","bibi signal",0.02, 0.0005, 0.1);
  RooRealVar Ts("Ts","tti signal",0.2,0.0005,1.);
  ms.setVal(1066.2);   ms.setConstant(kTRUE);
  ns.setVal(0.733);    ns.setConstant(kTRUE);
  n2s.setVal(0.95);   n2s.setConstant(kTRUE);
  bbs.setVal(0.0323); // bbs.setConstant(kTRUE);
  Ts.setVal(0.0535);   // Ts.setConstant(kTRUE);

  // RooTsallis* rt = new RooTsallis("rt","rt",pt,m,n,T);
  // ws->import(*rt);
  RooTsallis3* rt2 = new RooTsallis3("rt2","rt2",pt,ms,ns,n2s,bbs,Ts);
  ws->import(*rt2);
  RooTsallis3* rt3 = new RooTsallis3("rt3","rt3",pt,m,n,n2,bb,T);
  ws->import(*rt3);
  ws->factory("Gaussian::gau(pt,mean[3.,2.,15.],sigma[2.,0.1,10.])");
  ws->factory("Chebychev::che(pt,{a0[0.5,0.,1.],a1[0.5,0.,1.]})");
  // ws->factory("SUM::all(coeffPol[0.1,0.,1.]*che,coeffGau[0.1,0.,1.]*gau,rt)");
  ws->factory("SUM::all(coeffGau[0.1,0.,1.]*gau,rt2)");

  // signal fit
  RooFitResult* bkgfit = ws->pdf("rt3")->fitTo(*bkg,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  RooPlot *framebkg = pt.frame();

  char reducestr[300];
  sprintf(reducestr,"pt < %f && pt > %f",varMax,varMin);
  
  bkg->plotOn(framebkg,DataError(RooAbsData::SumW2),Cut(reducestr));
  ws->pdf("rt3")->plotOn(framebkg,LineColor(kBlue),Normalization(bkg->sumEntries(),RooAbsReal::NumEvent));
  RooHist *hpullbkg = framebkg->pullHist();
  RooHist *hresidbkg = framebkg->residHist();
 
  char fileToSave[200];
  sprintf(fileToSave,"text/paramsBkg_%dGeV_all.txt",int(mZZcenter));
  ofstream os1(fileToSave);
  (RooArgSet(bkgfit->floatParsFinal())).writeToStream(os1,false);
  os1.close();

  RooFitResult* sigfit = ws->pdf("rt2")->fitTo(*sig,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  RooPlot *framesig = pt.frame();

  sig->plotOn(framesig,DataError(RooAbsData::SumW2),Cut(reducestr));
  ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue),Normalization(sig->sumEntries(),RooAbsReal::NumEvent));
  RooHist *hpullsig = framesig->pullHist();
  RooHist *hresidsig = framesig->residHist();

  sprintf(fileToSave,"text/paramsSig_%dGeV_all.txt",int(mZZcenter));
  ofstream os2(fileToSave);
  (RooArgSet(sigfit->floatParsFinal())).writeToStream(os2,false);
  os2.close();

  TCanvas can("can","The canvas",5.,5.,1000.,1300.); 
  can.Divide(2,4);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);

  can.cd(1);
  gPad->SetBottomMargin(0.0);
  framesig->Draw();
  sprintf(fileToSave,"Signal %d GeV - HRes",int(mZZcenter));
  t->DrawLatex(0.6,0.8,fileToSave); 
  // t->DrawLatex(0.6,0.8,"Signal 125 GeV - HRes"); 
  can.cd(2);
  gPad->SetBottomMargin(0.0);
  framebkg->Draw();  
  sprintf(fileToSave,"ZZ around %d GeV - POWHEG",int(mZZcenter));
  t->DrawLatex(0.6,0.8,fileToSave); 

  can.cd(3);
  gPad->SetLogy();
  gPad->SetTopMargin(0.0);
  framesig->Draw();
  can.cd(4);
  gPad->SetLogy(); 
  gPad->SetTopMargin(0.0);
  framebkg->Draw();

  // Compute ratio hists
  double *xsig    = hresidsig->GetX();
  double *yressig = hresidsig->GetY();
  for (int i = 0; i < nbins2; i++) {
 
    pt.setVal(xsig[i]);
    float thisweight = sig->weight(RooArgSet(pt));
    double thiserrup, thiserrdown;
    sig->weightError(thiserrup, thiserrdown, RooAbsData::SumW2);
    // cout << "yres[" << i << "] = " << yressig[i] << " x = " << xsig[i] << " weight = " << thisweight << endl;
    yressig[i] = thisweight/(thisweight-yressig[i]);
    float yerrsigup = thiserrup/(thisweight-yressig[i]);
    float yerrsigdown = thiserrdown/(thisweight-yressig[i]);
    hresidsig->SetPoint(i,xsig[i],yressig[i]);
    hresidsig->SetPointError(i,0.,0.,yerrsigdown,yerrsigup);

  } 

  double *xbkg    = hresidbkg->GetX();
  double *yresbkg = hresidbkg->GetY();
  for (int i = 0; i < nbins; i++) {
 
    pt.setVal(xbkg[i]);
    float thisweight = bkg->weight(RooArgSet(pt));
    double thiserrup, thiserrdown;
    bkg->weightError(thiserrup, thiserrdown, RooAbsData::SumW2);
    yresbkg[i] = thisweight/(thisweight-yresbkg[i]);
    float yerrbkgup = thiserrup/(thisweight-yresbkg[i]);
    float yerrbkgdown = thiserrdown/(thisweight-yresbkg[i]);
    hresidbkg->SetPoint(i,xbkg[i],yresbkg[i]);
    hresidbkg->SetPointError(i,0.,0.,yerrbkgdown,yerrbkgup);

  }
 
  RooPlot* ressig = pt.frame(Title("Residuals Distribution")) ;
  ressig->GetYaxis()->SetTitle("Ratio");
  /* ressig->SetLabelSize(0.08,"XYZ");
  ressig->SetTitleSize(0.08,"XYZ");
  ressig->SetTitleOffset(0.6,"Y");
  ressig->SetTitleOffset(1.0,"X"); */
  ressig->addPlotable(hresidsig,"P") ; 
  ressig->SetMinimum(0.5); 
  ressig->SetMaximum(1.5); 

  can.cd(5);
  gPad->SetBottomMargin(0.0);
  gPad->SetGridy();
  ressig->Draw();

  RooPlot* resbkg = pt.frame(Title("Residuals Distribution")) ;
  resbkg->GetYaxis()->SetTitle("Ratio");
  /* resbkg->SetLabelSize(0.08,"XYZ");
  resbkg->SetTitleSize(0.08,"XYZ");
  resbkg->SetTitleOffset(0.6,"Y");
  resbkg->SetTitleOffset(1.0,"X"); */
  resbkg->addPlotable(hresidbkg,"P") ; 
  resbkg->SetMinimum(0.5); 
  resbkg->SetMaximum(1.5); 

  can.cd(6);
  gPad->SetBottomMargin(0.0);
  gPad->SetGridy();
  resbkg->Draw();

  RooPlot* pulsig = pt.frame(Title("Pull Distribution")) ;
  pulsig->GetYaxis()->SetTitle("Pull");
  /* pulsig->SetLabelSize(0.08,"XYZ");
  pulsig->SetTitleSize(0.08,"XYZ");
  pulsig->SetTitleOffset(0.6,"Y");
  pulsig->SetTitleOffset(1.0,"X"); */
  pulsig->addPlotable(hpullsig,"P") ; 
  pulsig->SetMinimum(-8.); 
  pulsig->SetMaximum(8.); 

  can.cd(7);
  gPad->SetTopMargin(0.0);
  gPad->SetGridy();
  pulsig->Draw();

  RooPlot* pulbkg = pt.frame(Title("Pull Distribution")) ;
  pulbkg->GetYaxis()->SetTitle("Pull");
  /* pulbkg->SetLabelSize(0.08,"XYZ");
  pulbkg->SetTitleSize(0.08,"XYZ");
  pulbkg->SetTitleOffset(0.6,"Y");
  pulbkg->SetTitleOffset(1.0,"X"); */
  pulbkg->addPlotable(hpullbkg,"P") ; 
  pulbkg->SetMinimum(-8.); 
  pulbkg->SetMaximum(8.); 

  can.cd(8);
  gPad->SetTopMargin(0.0);
  gPad->SetGridy();
  pulbkg->Draw();

  
  sprintf(fileToSave,"figs/fitSigBkg_%dGeV_all.pdf",int(mZZcenter));
  can.SaveAs(fileToSave);
 

}
