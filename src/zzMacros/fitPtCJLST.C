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

void Run2D(RooWorkspace *ws, int mZZcenter) {
  
  float minZZ = 105.;
  float maxZZ = 150.;

  char fileToSave[200];
  // TOYS TO ESTIMATE f_VBF UNCERTAINTY

  // Read data files
  TChain* sigTree = new TChain("SelectedTree");
  sigTree->Add("MELA/datafiles/4e/HZZ4lTree_H125.root");
  sigTree->Add("MELA/datafiles/4mu/HZZ4lTree_H125.root");
  sigTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_H125.root"); 

  TChain* bkgTree = new TChain("SelectedTree");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo2e2mu.root");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo2e2tau.root");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo2mu2tau.root");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo4e.root");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo4mu.root");
  bkgTree->Add("MELA/datafiles/4e/HZZ4lTree_ZZTo4tau.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo2e2mu.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo2e2tau.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo2mu2tau.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo4e.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo4mu.root");
  bkgTree->Add("MELA/datafiles/4mu/HZZ4lTree_ZZTo4tau.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo2e2mu.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo4e.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo4mu.root");
  bkgTree->Add("MELA/datafiles/2mu2e/HZZ4lTree_ZZTo4tau.root"); 

  /// Fit mass distributions
  RooRealVar ZZMass("ZZMass","ZZMass", 110.,minZZ,maxZZ,"GeV/c");
  ws->import(ZZMass);
  sprintf(fileToSave,"ZZMass > %f && ZZMass < %f",minZZ,maxZZ);
  RooDataSet* sigM = new RooDataSet("sigM","signal mass",sigTree,RooArgSet(ZZMass),fileToSave); 
  RooDataSet* bkgM = new RooDataSet("bkgM","ZZ mass",bkgTree,RooArgSet(ZZMass),fileToSave); 

  /* ws->factory("CBShape::sigCB1(ZZMass,meanSig[125.,105.,140.],sigmaSig1[2.,1.,10.],alpha[1.7,0.3,3.0],enne[2.,1.5,50.])");
  ws->factory("Gaussian::sigCB2(ZZMass,meanSig,sigmaSig2[8.,4.,15.])");
  ws->factory("SUM::sigCB(fCB1[0.5,0.001.,0.999]*sigCB1,sigCB2)"); */  
  RooDataHist* sigMH = new RooDataHist("sigMH","sigMH",RooArgSet(ZZMass),*sigM); 
  RooHistPdf* sigCB = new RooHistPdf("sigCB","sigCB",RooArgSet(ZZMass),*sigMH);
  ws->import(*sigCB);
  ws->factory("Chebychev::bkgCh(ZZMass,{a0[0.1,-1.,1.],a1[0.1,-1.,1.],a2[0.1,-1.,1.]})");

  // ws->pdf("sigCB")->fitTo(*sigM,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));
  ws->pdf("bkgCh")->fitTo(*bkgM,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1)); 

  RooPlot *framebkgM = ZZMass.frame();

  bkgM->plotOn(framebkgM,DataError(RooAbsData::SumW2));
  ws->pdf("bkgCh")->plotOn(framebkgM,LineColor(kBlue),Normalization(bkgM->sumEntries(),RooAbsReal::NumEvent));
  
  RooPlot *framesigM = ZZMass.frame();
 
  sigM->plotOn(framesigM,DataError(RooAbsData::SumW2));
  ws->pdf("sigCB")->plotOn(framesigM,LineColor(kBlue),Normalization(sigM->sumEntries(),RooAbsReal::NumEvent));
  
  TCanvas can2("can2","The canvas",5.,5.,800.,500.); 
  can2.Divide(2,1);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);

  can2.cd(1);  framesigM->Draw();
  can2.cd(2);  framebkgM->Draw();

  sprintf(fileToSave,"figs/fitMass_%dGeV.pdf",int(mZZcenter));
  can2.SaveAs(fileToSave);

  /// Fix all shapes
  ws->var("a0")->setConstant(kTRUE);
  ws->var("a1")->setConstant(kTRUE);
  ws->var("a2")->setConstant(kTRUE);
  ws->var("m")->setConstant(kTRUE);
  ws->var("n")->setConstant(kTRUE);
  ws->var("bb")->setConstant(kTRUE);
  ws->var("n2")->setConstant(kTRUE);
  ws->var("bb2")->setConstant(kTRUE);
  ws->var("T")->setConstant(kTRUE);
  ws->var("fexp")->setConstant(kTRUE);
  ws->var("ms")->setConstant(kTRUE);
  ws->var("ns")->setConstant(kTRUE);
  ws->var("bbs")->setConstant(kTRUE);
  ws->var("n2s")->setConstant(kTRUE);
  ws->var("bb2s")->setConstant(kTRUE);
  ws->var("Ts")->setConstant(kTRUE);
  ws->var("fexps")->setConstant(kTRUE);
  ws->var("mv")->setConstant(kTRUE);
  ws->var("nv")->setConstant(kTRUE);
  ws->var("bbv")->setConstant(kTRUE);
  ws->var("n2v")->setConstant(kTRUE);
  ws->var("bb2v")->setConstant(kTRUE);
  ws->var("Tv")->setConstant(kTRUE);
  ws->var("fexpv")->setConstant(kTRUE);

  // Generate and fit toys

  //1D
  ws->factory("SUM::ptSig(fVBF[0.5,0.,1.]*rt,rt2)");
  ws->factory("SUM::massPdf(fBkg[0.5,0.,1.]*bkgCh,sigCB)");
  ws->factory("SUM::ptPdf(fBkg*rt4,ptSig)");

  //2D
  ws->factory("PROD::allGG(sigCB,rt2)");  
  ws->factory("PROD::allVBF(sigCB,rt)");
  ws->factory("PROD::allZZ(bkgCh,rt4)");  
  ws->factory("PROD::allSig(sigCB,ptSig)");
  ws->factory("SUM::all(fBkg*allZZ,allSig)");

  ifstream theFile("configToys.txt");
  char thePar[10];
  float theVal;	
  int nToys;
  float fBkgVal, fVBFVal, lumiVal;

  // fBkg and fVBF from Mike's files for 20 fb-1: can be overwritten
  // *** ee
  // process	ttH	ZH	WH	qqH	ggH	QQZZ	GGZZ	FAKES	
  // rate	0.0109615892	0.0400934661727	0.0727344736236	0.17448742481	 1.94373327139	3.44729957187	0.0691200171852	0.966791503492	
  // *** emu
  // process	ttH	ZH	WH	qqH	ggH	QQZZ	GGZZ	FAKES	
  // rate	0.0297855894832	0.108944743558	0.197639150046	0.449271834373	 5.28164668738	9.12046958147	0.171200846394	2.409869007	
  // ** mumu
  // process	ttH	ZH	WH	qqH	ggH	QQZZ	GGZZ	FAKES	
  // rate	0.0208474338004	0.0762522538126	0.138330956964	0.334468868705	 3.69671312815	7.94918247279	0.137969786197	1.51563994758	
  
  fVBFVal = (0.1745+0.4493+0.3345)/(0.1745+0.4493+0.3345+2.0667+5.6201+3.9321);  
  float allEvents = 4.4822+11.7015+9.6028+0.1745+0.4493+0.3345+2.0667+5.6201+3.9321;  
  fBkgVal = (4.4822+11.7015+9.6028)/allEvents;
  float scaleBkg = (maxZZ-minZZ)/35.;  // Mike's numbers are for 105-140 GeV
  fBkgVal = scaleBkg*fBkgVal/(1. + fBkgVal*(scaleBkg-1.));

  while (theFile >> thePar >> theVal) {
    if (!strcmp(thePar,"n")) {nToys = int(theVal);}
    if (!strcmp(thePar,"fVBF")) {fVBFVal = theVal;}
    if (!strcmp(thePar,"fBkg")) {fBkgVal = theVal;}
    if (!strcmp(thePar,"lumi")) {lumiVal = theVal;}
  }

  cout << "fBkg GEN = " << fBkgVal << endl;
  cout << "fVBF GEN = " << fVBFVal << endl;
  cout << "lumi = " << lumiVal << " fb-1" << endl;

  float howManyEvents = allEvents*lumiVal/20.;
  cout << "Total number of expected events = " << howManyEvents << endl;
  cout << "Number of expected signal events = " << howManyEvents*(1-fBkgVal) << endl;
  int howManyEventsGG = int(howManyEvents*(1-fBkgVal)*(1-fVBFVal));
  int howManyEventsVBF = int(howManyEvents*(1-fBkgVal)*fVBFVal);
  int howManyEventsBKG = int(howManyEvents*fBkgVal);
  
  // Store results
  TH1F* fVBFRes = new TH1F("fVBFRes","fVBF",11,-0.05,1.05);
  TH1F* fVBFErr = new TH1F("fVBFErr","fVBF",10,-0.05,0.4);
  TH1F* fVBFPull = new TH1F("fVBFPull","fVBF",10,-5.,5.);
  TH1F* nll = new TH1F("nll","fVBF",11,-100000.,100000.);

  for (unsigned int iToy = 0; iToy < nToys; iToy++) {

    cout << endl << "####" << endl;
    cout << "Generating toy experiment n. " << iToy+1 << endl;

    ws->var("fVBF")->setVal(fVBFVal);
    ws->var("fVBF")->setConstant(kFALSE);
    ws->var("fBkg")->setVal(fBkgVal);
    ws->var("fBkg")->setConstant(kFALSE);

    TDatime *now = new TDatime();
    Int_t seed = now->GetDate() + now->GetTime();
    cout << "RooFit Generation Seed = " << seed << endl;
    RooRandom::randomGenerator()->SetSeed(seed);
    cout << "####" << endl << endl;

    RooDataSet *dataToy = ws->pdf("allGG")->generate(RooArgSet(*(ws->var("ZZMass")),*(ws->var("pts"))),howManyEventsGG,Extended());
    RooDataSet *dataToy2 = ws->pdf("allVBF")->generate(RooArgSet(*(ws->var("ZZMass")),*(ws->var("pts"))),howManyEventsVBF,Extended());
    RooDataSet *dataToy3 = ws->pdf("allZZ")->generate(RooArgSet(*(ws->var("ZZMass")),*(ws->var("pts"))),howManyEventsBKG,Extended());
    dataToy->append(*dataToy2);
    dataToy->append(*dataToy3);

    // int nFitPar;  
    Double_t theNLL;
    ws->pdf("massPdf")->fitTo(*dataToy,Minos(0),SumW2Error(kTRUE),NumCPU(1));  
    ws->var("fBkg")->setConstant(kTRUE);
  
    // RooFitResult *rfr = ws->pdf("all")->fitTo(*dataToy,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));      
    RooFitResult *rfr = ws->pdf("ptPdf")->fitTo(*dataToy,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));

    // nFitPar = rfr->floatParsFinal().getSize();  
    theNLL = rfr->minNll();  
    float fVBF_static = ws->var("fVBF")->getVal();
    float fVBFErr_static = ws->var("fVBF")->getError();
    fVBFRes->Fill(fVBF_static);
    fVBFErr->Fill(fVBFErr_static);
    fVBFPull->Fill(fVBF_static/fVBFErr_static);
    nll->Fill(theNLL);

    if (iToy == 0) { // JUST DRAW THE FIRST EXPERIMENT

      RooHist *hresid = new RooHist();
      cout << "OK qui" << endl;

      RooPlot *mframe = ws->var("ZZMass")->frame();

      TCanvas can3("can3","The canvas",5.,5.,1000.,800.); 
      can3.Divide(2,2);

      dataToy->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(14));      
      ws->pdf("all")->plotOn(mframe,Components("allZZ"),LineColor(kBlue),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      ws->pdf("all")->plotOn(mframe,Components("allZZ,allVBF"),LineColor(kRed),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      ws->pdf("all")->plotOn(mframe,LineColor(kBlack),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      hresid = mframe->pullHist();

      RooPlot* mframeres = ws->var("ZZMass")->frame(Title("Residuals Distribution")) ;
      mframeres->addPlotable(hresid,"P") ;  
      
      can3.cd(1);mframe->Draw();
      can3.cd(3);mframeres->Draw();
      
      RooBinning myRb(ws->var("pts")->getMin(), ws->var("pts")->getMax());
      myRb.addBoundary(2.);
      myRb.addBoundary(5.);
      myRb.addBoundary(7.5);
      myRb.addBoundary(10.);
      myRb.addBoundary(12.5);
      myRb.addBoundary(15.);
      myRb.addBoundary(18.);
      myRb.addBoundary(30.);
      myRb.addBoundary(50.);
      myRb.addBoundary(100.);
      myRb.addBoundary(200.);
    
      RooPlot *tframe = ws->var("pts")->frame();      
      
      dataToy->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(myRb));
      ws->pdf("all")->plotOn(tframe,Components("allZZ"),LineColor(kBlue),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      ws->pdf("all")->plotOn(tframe,Components("allZZ,allVBF"),LineColor(kRed),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      ws->pdf("all")->plotOn(tframe,LineColor(kBlack),Normalization(dataToy->sumEntries(),RooAbsReal::NumEvent));
      hresid = tframe->pullHist();

      RooPlot* tframeres =  ws->var("pts")->frame(Title("Residuals Distribution")) ;
      tframeres->addPlotable(hresid,"P") ;  
	
      can3.cd(2);gPad->SetLogx();gPad->SetLogy();tframe->Draw();
      can3.cd(4);gPad->SetLogx();tframeres->Draw();

      sprintf(fileToSave,"figs/fit2D_%dGeV.pdf",int(mZZcenter));
      can3.SaveAs(fileToSave);

    } 

  }

  TCanvas can4("can4","The canvas",5.,5.,1000.,800.); 
  can4.Divide(2,2);
  can4.cd(1);fVBFRes->Draw();
  can4.cd(2);fVBFErr->Draw();
  can4.cd(3);fVBFPull->Draw();
  can4.cd(4);nll->Draw();
  sprintf(fileToSave,"figs/toyResults_%dGeV.pdf",int(mZZcenter));
  can4.SaveAs(fileToSave);

  return; 
}

void fitPtCJLST(int LHCsqrts = 7, int whichtype = 1)

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
  TH1F* ptH = (TH1F*)input.Get("ptH");
  if (rebinType[whichtype] > 1) ptH->Rebin(rebinType[whichtype]);
  if (maxType[whichtype] < ptH->GetBinLowEdge(ptH->GetNbinsX() + 1) - ptH->GetBinWidth(1)) {
    int theBin = ptH->FindBin(maxType[whichtype]);
    ptH->GetXaxis()->SetRange(1,theBin-1);
  }

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // cout << endl << "Signal " << endl;   

  RooDataHist* rdh = new RooDataHist("rdh","Some dataset",RooArgList(*pt),Import(*ptH,kFALSE));
 
  // fit definitions
  RooWorkspace *ws = new RooWorkspace("ws");

  RooRealVar m("m","emme", 110.,10., 600.,"GeV/c");
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
      n2.setVal(1.0678);   n2.setConstant(kTRUE);
      n.setVal(1.010);   n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.20);   // T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    } else {
      m.setVal(83.54);   // m.setConstant(kTRUE);    
      bb.setVal(0.020); // bb.setConstant(kTRUE);
      n2.setVal(1.0803);   n2.setConstant(kTRUE);
      n.setVal(1.0207);    n.setConstant(kTRUE);
      bb2.setVal(100000.);  bb2.setConstant(kTRUE);
      T.setVal(0.20);   // T.setConstant(kTRUE);
      fexp.setVal(0.0);    fexp.setConstant(kTRUE);
    }
  }
  else if (whichtype == 1) {
    m.setVal(235.3);   // m.setConstant(kTRUE);
    n.setVal(1.1705);   n.setConstant(kTRUE);
    n2.setVal(4.086);   n2.setConstant(kTRUE);
    bb.setVal(0.0294); // bb.setConstant(kTRUE);
    T.setVal(0.000109);   // T.setConstant(kTRUE);
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
      n.setVal(0.6048);   n.setConstant(kTRUE);
      n2.setVal(0.9086);   n2.setConstant(kTRUE);
      bb.setVal(0.0280);   // bb.setConstant(kTRUE);
      T.setVal(0.0866);   // T.setConstant(kTRUE);
      bb2.setVal(0.00657);   bb2.setConstant(kTRUE);
      fexp.setVal(0.0849);   fexp.setConstant(kTRUE);
    } else {
      m.setVal(456.99);   // m.setConstant(kTRUE);
      n.setVal(1.068);   n.setConstant(kTRUE);
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

  char fileToSave[200];
  sprintf(fileToSave,"text/paramsCJLST_%s_%dTeV_all.txt",nameSample[whichtype].c_str(),LHCsqrts);
  ofstream os1(fileToSave);
  (RooArgSet(fit->floatParsFinal())).writeToStream(os1,false);
  (RooArgSet(fit->constPars())).writeToStream(os1,false);
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

  sprintf(fileToSave,"figs/fitCJLST_%s_%dTeV.pdf",nameSample[whichtype].c_str(),LHCsqrts);
  can.SaveAs(fileToSave);

}



