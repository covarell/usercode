//Fit the PT curves for sig and BKG
//
//By: Roberto Covarelli (U of Rochester)

#include <math.h>
#include <string>
#include <vector>
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
#include "RooBinning.h"
#include "RooFitResult.h"

// #include "PDFs/RooTsallis.h"
// #include "PDFs/RooTsallis2.h"
#include "PDFs/RooTsallis3.h" 

using namespace RooFit;

RooDataHist* withSmartBinning(TH1F* source, RooRealVar* var, float min, float max, int type) {
  
  // DEFINE BINNING
  // type:
  // 0 = background
  // 1 = gg signal
  // 2 = VBF signal
  float lowlimit;
  std::vector<float> sizes;
  std::vector<float> highlimits;  
  if (type == 0) {
    lowlimit = min;
    sizes.push_back(1.);
    highlimits.push_back(max);
  } else if (type == 1) {
    lowlimit = min;
    sizes.push_back(4.);        
    // sizes.push_back(8.);        
    // sizes.push_back(4.);
    // sizes.push_back(8.);
    // highlimits.push_back(150.); 
    // highlimits.push_back(200.);
    // highlimits.push_back(320.);
    highlimits.push_back(max);
  } else {
    lowlimit = min;
    sizes.push_back(4.);
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
  float binWidth = source->GetBinWidth(1);

  float thisWeight = 0.;
  float thisWeightsum = 0.;
  int j = 0;
  int whichInterval = 0;
  for (Int_t i=1; i<=source->GetNbinsX(); i++) {

    float thispt = source->GetXaxis()->GetBinCenter(i);
    if (thispt > lowlimit && thispt < highlimits.back()) {
      if (thispt > highlimits.at(whichInterval)) whichInterval++;
      var->setVal(thispt);
      thisWeight += source->GetBinContent(i)*nWeight;
      thisWeightsum += source->GetBinError(i)*nWeight*source->GetBinError(i)*nWeight;
      j++;
      cout << "Entries in TH1 bin " << i << " = " << source->GetBinContent(i)*nWeight << " +/- " << source->GetBinError(i)*nWeight << endl;   
      if (binWidth*j >= sizes.at(whichInterval)) {
	cout << "Entro" << endl; 
	if (thisWeight <= 0) thisWeight = 0.1; 
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

TH1F* NNLOLowPtReweight(TH1F* toReweight, TH1F* weighter, float upTo) {
  // normalize so that weight is 1 at the maximum reweighting point
  int binUpTo = toReweight->FindBin(upTo);
  float scale = toReweight->GetBinContent(binUpTo)/weighter->GetBinContent(binUpTo);
  weighter->Scale(scale);
  for (int i = 1; i <=binUpTo; i++) {  // low pt
    toReweight->SetBinContent(i,weighter->GetBinContent(i));
    toReweight->SetBinError(i,sqrt(weighter->GetBinError(i)*weighter->GetBinError(i) + weighter->GetBinError(i)*weighter->GetBinError(i)));
  }
  return toReweight; 
}

void fitPt(float mZZcenter = 126., float mZZspread = 5., int LHCsqrts = 7, bool isVBFsignal = false)
{

  float varMinBkg = 1.;        // all in GeV
  float varMaxBkg = 180.;
  float varMinSig = 3.;
  float varMaxSig = 400.;

  char fileToOpen[200];
  sprintf(fileToOpen,"PT_Y_%dTeV.root",LHCsqrts);
  TFile* fileb = new TFile(fileToOpen);
  // sprintf(fileToOpen,"PT_Y_gg125-200_%dTeV.root",LHCsqrts);
  sprintf(fileToOpen,"PT_Y_%dTeV.root",LHCsqrts);
  TFile* files = new TFile(fileToOpen);
  
  TH1F* massH = (TH1F*)((TH2F*)fileb->Get("Pt_bkg"))->ProjectionX();
  float mZZmin = mZZcenter - mZZspread;
  float mZZmax = mZZcenter + mZZspread;
  int binMin = massH->FindBin(mZZmin);
  int binMax = massH->FindBin(mZZmax);

  // Check strange way of filling histos (CM)
  int binMean = (binMin+binMax)/2;
  if ((fabs(massH->GetBinContent(binMean) - 1.) < 0.01 || 
       massH->GetBinContent(binMean) == 0.) &&
      (fabs(massH->GetBinContent(binMean-1) - 1.) < 0.01 || 
       massH->GetBinContent(binMean-1) == 0.) &&
      (fabs(massH->GetBinContent(binMean+1) - 1.) < 0.01 || 
       massH->GetBinContent(binMean+1) == 0.)
      ) {
    binMin = binMean;    binMax = binMean;
  }

  TH1F* bkgH = (TH1F*)((TH2F*)fileb->Get("Pt_bkg"))->ProjectionY("bkgH",binMin,binMax);
  TH1F* sigH;
  
  if (isVBFsignal) {
    sprintf(fileToOpen,"PT_Y_VBF_%dTeV.root",LHCsqrts); 
    TFile* files2 = new TFile(fileToOpen);
    sigH = (TH1F*)((TH2F*)files2->Get("Pt_sigVBF"))->ProjectionY("sigH",binMin,binMax); 
  } else sigH = (TH1F*)((TH2F*)files->Get("Pt_sig"))->ProjectionY("sigH",binMin,binMax);
 
  /* sprintf(fileToOpen,"PT_%d_Temp.root",int(mZZcenter));
  TFile* file2 = new TFile(fileToOpen);
  TH1F* sigH = (TH1F*)file2->Get("Pt_sig");  */
   
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // Build datasets
  RooRealVar* pt = new RooRealVar("pt","p_{T}^{H}",varMinBkg,varMaxBkg,"GeV/c");
  
  // TH1F* bkgHtest = new TH1F("bkgHtest","bkgHtest",nbins,varMinBkg,varMaxBkg);
  // RooDataHist* bkg = new RooDataHist("bkg","Background dataset",RooArgList(pt),Import(*bkgH,kTRUE),Weight(nBkgWeight));
  
  cout << endl << "Background " << endl; 
  RooDataHist* bkg = withSmartBinning(bkgH,pt,varMinBkg,varMaxBkg,0);
  bkg->SetName("bkg");

  for (Int_t i=0; i<bkg->numEntries(); i++) {

    const RooArgSet* aRow = bkg->get(i);
    RooRealVar* ptprime = (RooRealVar*)aRow->find("pt");
    pt->setVal(ptprime->getVal());
    cout << "Entries in RooDataSet bin " << i << " = " << bkg->weight() << " +/- " << bkg->weightError(RooAbsData::SumW2) << endl;  
  } 

  RooRealVar* pts = new RooRealVar("pts","p_{T}^{H}",varMinSig,varMaxSig,"GeV/c");

  cout << endl << "Signal " << endl;   

  // RooDataHist* sig = new RooDataHist("sig","Signal dataset",RooArgList(pts),Import(*sigH,kTRUE),Weight(nSigWeight));

  int whichSig = 1;
  if (isVBFsignal) whichSig = 2;
  RooDataHist* sig = withSmartBinning(sigH,pts,varMinSig,varMaxSig,whichSig);
  sig->SetName("sig");
 
  for (Int_t i=0; i<sig->numEntries(); i++) {

    const RooArgSet* aRow = sig->get(i);
    RooRealVar* ptprime = (RooRealVar*)aRow->find("pts");
    pts->setVal(ptprime->getVal());
    cout << "Entries in RooDataSet bin " << i << " = " << sig->weight() << " +/- " << sig->weightError(RooAbsData::SumW2) << endl;  
  } 

  // fit definitions
  RooWorkspace *ws = new RooWorkspace("ws");

  RooRealVar m("m","emme", 110.,10., 600.,"GeV/c");
  RooRealVar n("n","enne", 0.93, 0.5, 15.);
  RooRealVar n2("n2","enne2", 0.75, 0.5, 15.);
  RooRealVar bb("bb","bibi",0.02, 0.0005, 0.1);
  RooRealVar T("T","tti",0.2,0.0005,1.);
  RooRealVar bb2("bb2","bibi2",0.02, 0.0005, 0.1);
  RooRealVar fexp("fexp","f_exp",0.02, 0.0, 1.0);
  m.setVal(54.47);   m.setConstant(kTRUE);
  n.setVal(1.255);   // n.setConstant(kTRUE);
  n2.setVal(1.73);   n2.setConstant(kTRUE);
  bb.setVal(0.020); // bb.setConstant(kTRUE);
  bb2.setVal(0.020);  bb2.setConstant(kTRUE);
  T.setVal(0.20);   // T.setConstant(kTRUE);
  fexp.setVal(0.0);    fexp.setConstant(kTRUE);

  RooRealVar ms("ms","emme signal", 110.,10., 6000.,"GeV/c");
  RooRealVar ns("ns","enne signal", 0.93, 0.5, 15.);
  RooRealVar n2s("n2s","enne2 signal", 0.75, 0.5, 15.); 
  RooRealVar bbs("bbs","bibi signal",0.02, 0.0005, 0.1);
  RooRealVar Ts("Ts","tti signal",0.2,0.00000005,1.);
  RooRealVar bb2s("bb2s","bibi2 signal",0.02, 0.0005, 0.1);
  RooRealVar fexps("fexps","f_exp signal",0.02, 0.0, 1.0);
  if (isVBFsignal) {
    ms.setVal(653.8);   ms.setConstant(kTRUE);
    ns.setVal(1.1705);   ns.setConstant(kTRUE);
    n2s.setVal(4.086);   n2s.setConstant(kTRUE);
    bbs.setVal(0.0294); // bbs.setConstant(kTRUE);
    Ts.setVal(0.0064);   // Ts.setConstant(kTRUE);
    bb2s.setVal(0.020);   // bb2.setConstant(kTRUE);
    fexps.setVal(0.1);   // fexp.setConstant(kTRUE);
  } else {
    ms.setVal(1803.8);   ms.setConstant(kTRUE);
    ns.setVal(0.733);   ns.setConstant(kTRUE);
    n2s.setVal(0.95);   n2s.setConstant(kTRUE);
    bbs.setVal(0.0323); // bbs.setConstant(kTRUE);
    Ts.setVal(0.064);   // Ts.setConstant(kTRUE);
    bb2s.setVal(0.020);   // bb2.setConstant(kTRUE);
    fexps.setVal(0.1);   // fexp.setConstant(kTRUE);
  }
  
  // RooTsallis* rt = new RooTsallis("rt","rt",pt,m,n,T);
  // ws->import(*rt);
  RooTsallis3* rt2 = new RooTsallis3("rt2","rt2",*pts,ms,ns,n2s,bbs,bb2s,Ts,fexps);
  ws->import(*rt2);
  RooTsallis3* rt3 = new RooTsallis3("rt3","rt3",*pt,m,n,n2,bb,bb2,T,fexp);
  ws->import(*rt3);
  ws->factory("Gaussian::gau(pt,mean[3.,2.,15.],sigma[2.,0.1,10.])");
  ws->factory("Chebychev::che(pt,{a0[0.5,0.,1.],a1[0.5,0.,1.]})");
  // ws->factory("SUM::all(coeffPol[0.1,0.,1.]*che,coeffGau[0.1,0.,1.]*gau,rt)");
  ws->factory("SUM::all(coeffGau[0.1,0.,1.]*gau,rt2)");

  // signal fit
  RooFitResult* bkgfit = ws->pdf("rt3")->fitTo(*bkg,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  RooPlot *framebkg = pt->frame();

  char reducestr[300];
  sprintf(reducestr,"pt > %f && pt < %f",pt->getMin(),pt->getMax());
  
  bkg->plotOn(framebkg,DataError(RooAbsData::SumW2),Cut(reducestr));
  ws->pdf("rt3")->plotOn(framebkg,LineColor(kBlue),Normalization(bkg->sumEntries(),RooAbsReal::NumEvent));
  RooHist *hpullbkg = framebkg->pullHist();
  RooHist *hresidbkg = framebkg->residHist();
 
  char fileToSave[200];
  sprintf(fileToSave,"text/paramsBkg_%dGeV_%dTeV_all.txt",int(mZZcenter),LHCsqrts);
  ofstream os1(fileToSave);
  (RooArgSet(bkgfit->floatParsFinal())).writeToStream(os1,false);
  os1.close();

  RooFitResult* sigfit = ws->pdf("rt2")->fitTo(*sig,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));  

  RooPlot *framesig = pts->frame();
  sprintf(reducestr,"pts > %f && pts < %f",pts->getMin(),pts->getMax());

  sig->plotOn(framesig,DataError(RooAbsData::SumW2),Cut(reducestr));
  // ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue),Normalization(sig->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue));
  RooHist *hpullsig = framesig->pullHist();
  RooHist *hresidsig = framesig->residHist();

  if (isVBFsignal) sprintf(fileToSave,"text/paramsVbf_%dGeV_%dTeV_all.txt",int(mZZcenter),LHCsqrts); 
  else sprintf(fileToSave,"text/paramsSig_%dGeV_%dTeV_all.txt",int(mZZcenter),LHCsqrts);

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
  if (isVBFsignal) sprintf(fileToSave,"VBF signal %d GeV - POWHEG",int(mZZcenter));
  else sprintf(fileToSave,"Signal %d GeV - HRes",int(mZZcenter));
  t->DrawLatex(0.6,0.8,fileToSave); 
  // t->DrawLatex(0.6,0.8,"Signal 125 GeV - HRes"); 
  can.cd(2);
  gPad->SetBottomMargin(0.0);
  framebkg->Draw();
  // gPad->SetLogy(); 
  // bkgHtest->Draw();
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
  int nbins2 = pts->getBinning().numBins();
  for (int i = 0; i < nbins2; i++) {
 
    pts->setVal(xsig[i]);
    float thisweight = sig->weight(RooArgSet(*pts));
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
  int nbins = pt->getBinning().numBins();
  for (int i = 0; i < nbins; i++) {
 
    pt->setVal(xbkg[i]);
    float thisweight = bkg->weight(RooArgSet(*pt));
    double thiserrup, thiserrdown;
    bkg->weightError(thiserrup, thiserrdown, RooAbsData::SumW2);
    yresbkg[i] = thisweight/(thisweight-yresbkg[i]);
    float yerrbkgup = thiserrup/(thisweight-yresbkg[i]);
    float yerrbkgdown = thiserrdown/(thisweight-yresbkg[i]);
    hresidbkg->SetPoint(i,xbkg[i],yresbkg[i]);
    hresidbkg->SetPointError(i,0.,0.,yerrbkgdown,yerrbkgup);

  }
 
  RooPlot* ressig = pts->frame(Title("Residuals Distribution")) ;
  ressig->GetYaxis()->SetTitle("Ratio");
  /* ressig->SetLabelSize(0.08,"XYZ");
  ressig->SetTitleSize(0.08,"XYZ");
  ressig->SetTitleOffset(0.6,"Y");
  ressig->SetTitleOffset(1.0,"X"); */
  ressig->addPlotable(hresidsig,"P") ; 
  ressig->SetMinimum(-0.5); 
  ressig->SetMaximum(2.5); 

  can.cd(5);
  gPad->SetBottomMargin(0.0);
  gPad->SetGridy();
  ressig->Draw();

  RooPlot* resbkg = pt->frame(Title("Residuals Distribution")) ;
  resbkg->GetYaxis()->SetTitle("Ratio");
  /* resbkg->SetLabelSize(0.08,"XYZ");
  resbkg->SetTitleSize(0.08,"XYZ");
  resbkg->SetTitleOffset(0.6,"Y");
  resbkg->SetTitleOffset(1.0,"X"); */
  resbkg->addPlotable(hresidbkg,"P") ; 
  resbkg->SetMinimum(-0.5); 
  resbkg->SetMaximum(2.5); 

  can.cd(6);
  gPad->SetBottomMargin(0.0);
  gPad->SetGridy();
  resbkg->Draw();

  RooPlot* pulsig = pts->frame(Title("Pull Distribution")) ;
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

  RooPlot* pulbkg = pt->frame(Title("Pull Distribution")) ;
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

  
  if (isVBFsignal) sprintf(fileToSave,"figs/fitVbfBkg_%dGeV_%dTeV_all.pdf",int(mZZcenter),LHCsqrts);
  else sprintf(fileToSave,"figs/fitSigBkg_%dGeV_%dTeV_all.pdf",int(mZZcenter),LHCsqrts);
  can.SaveAs(fileToSave);
 

}
