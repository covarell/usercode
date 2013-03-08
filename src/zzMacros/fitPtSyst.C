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

void fitPtSyst(float mZZcenter = 126., float mZZspread = 5., 
	       int LHCsqrts = 7, 
	       bool isVBFsignal = false, 
	       bool writeWeightHisto = false, 
	       bool run2D = false, string typeSyst = "Default")
{

  float varMinBkg = 1.;        // all in GeV
  float varMaxBkg = 142. + (mZZcenter/5.);
  if (LHCsqrts == 8 && fabs(mZZcenter - 125.) < 0.1) varMaxBkg = 170.;
  float varMinSig = 1.;
  float varMaxSig = 401.;
  float sizeSig = 4.;
  float minWeight = 0.1;

  char fileToOpen[200];
  sprintf(fileToOpen,"PT_Y_%dTeV.root",LHCsqrts);
  TFile* fileb = new TFile(fileToOpen);
  // sprintf(fileToOpen,"PT_Y_%dTeV.root",LHCsqrts);
  // sprintf(fileToOpen,"HRes/HResSystemtics_125.root");
  sprintf(fileToOpen,"HRes/HResSystv2_125.root");
  if (mZZcenter > 190.) sprintf(fileToOpen,"HRes/HResSystv2_200.root");
  if (mZZcenter > 390.) sprintf(fileToOpen,"HRes/HResSystv2_400.root");
  TFile* files = new TFile(fileToOpen);
  TFile* filegg;
  if (writeWeightHisto) {
    // sprintf(fileToOpen,"PT_Y_gg125-200_%dTeV.root",LHCsqrts);
    sprintf(fileToOpen,"PT_Y_gg125withPythia_%dTeV.root",LHCsqrts);
    if (mZZcenter > 190.) sprintf(fileToOpen,"PT_Y_gg200withPythia_%dTeV.root",LHCsqrts);
    if (mZZcenter > 390.) sprintf(fileToOpen,"PT_Y_gg400withPythia_%dTeV.root",LHCsqrts);
    filegg = new TFile(fileToOpen);  
  }

  TH1F* massH = (TH1F*)((TH2F*)fileb->Get("Pt_bkg"))->ProjectionX();
  float mZZmin = mZZcenter - mZZspread;
  float mZZmax = mZZcenter + mZZspread;
  int binMin = massH->FindBin(mZZmin);
  int binMax = massH->FindBin(mZZmax);

  RooRealVar* pt = new RooRealVar("pt","p_{T}^{H}",varMinBkg,varMaxBkg,"GeV/c");
  RooRealVar* pts = new RooRealVar("pts","p_{T}^{H}",varMinSig,varMaxSig,"GeV/c");
 
  TH1F* ggH = new TH1F();
  RooDataHist* gg = new RooDataHist();
  if (writeWeightHisto) {
    ggH = (TH1F*)((TH2F*)filegg->Get("Pt_sig"))->ProjectionY("ggH",binMin,binMax);
    gg = withSmartBinning(ggH,pts,varMinSig,varMaxSig,minWeight,sizeSig,2);
    gg->SetName("gg");
  }
  
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
  TH1F* sigVBFH;
  
  if (isVBFsignal) {
    sprintf(fileToOpen,"PT_Y_VBF_%dTeV.root",LHCsqrts); 
    TFile* files2 = new TFile(fileToOpen);
    sigVBFH = (TH1F*)((TH2F*)files2->Get("Pt_sigVBF"))->ProjectionY("sigVBFH",binMin,binMax); 
  } 
  // sigH = (TH1F*)((TH2F*)files->Get("Pt_sig"))->ProjectionY("sigH",binMin,binMax);
  sprintf(fileToOpen,"pt_%s",typeSyst.c_str());
  sigH = (TH1F*)(files->Get(fileToOpen));
 
  /* sprintf(fileToOpen,"PT_%d_Temp.root",int(mZZcenter));
  TFile* file2 = new TFile(fileToOpen);
  TH1F* sigH = (TH1F*)file2->Get("Pt_sig");  */
   
  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  
  // Build datasets
  
  // TH1F* bkgHtest = new TH1F("bkgHtest","bkgHtest",nbins,varMinBkg,varMaxBkg);
  // RooDataHist* bkg = new RooDataHist("bkg","Background dataset",RooArgList(pt),Import(*bkgH,kTRUE),Weight(nBkgWeight));
  
  cout << endl << "Background " << endl; 
  RooDataHist* bkg = withSmartBinning(bkgH,pt,varMinBkg,varMaxBkg,minWeight,sizeSig,0);
  bkg->SetName("bkg");

  for (Int_t i=0; i<bkg->numEntries(); i++) {

    const RooArgSet* aRow = bkg->get(i);
    RooRealVar* ptprime = (RooRealVar*)aRow->find("pt");
    pt->setVal(ptprime->getVal());
    cout << "Entries in RooDataSet bin " << i << " = " << bkg->weight() << " +/- " << bkg->weightError(RooAbsData::SumW2) << endl;  
  } 

 
  cout << endl << "Signal " << endl;   

  // RooDataHist* sig = new RooDataHist("sig","Signal dataset",RooArgList(pts),Import(*sigH,kTRUE),Weight(nSigWeight));
  RooDataHist* sigVBF = new RooDataHist();
  if (isVBFsignal) {
    sigVBF = withSmartBinning(sigVBFH,pts,varMinSig,varMaxSig,minWeight,sizeSig,3);
    sigVBF->SetName("sigVBF");
  }

  if (LHCsqrts == 7 && mZZcenter > 299. && mZZcenter < 301.) {
    sizeSig = 16.;
    minWeight = 0.005;
  }
  RooDataHist* sig = withSmartBinning(sigH,pts,varMinSig,varMaxSig,minWeight,sizeSig,1);
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
  bb2.setVal(0.20);  bb2.setConstant(kTRUE);
  T.setVal(0.20);   // T.setConstant(kTRUE);
  fexp.setVal(0.0);    fexp.setConstant(kTRUE);

  RooRealVar ms("ms","emme signal", 110.,10., 6000.,"GeV/c");
  RooRealVar ns("ns","enne signal", 0.93, 0.5, 15.);
  RooRealVar n2s("n2s","enne2 signal", 0.75, 0.5, 15.); 
  RooRealVar bbs("bbs","bibi signal",0.02, 0.0005, 0.1);
  RooRealVar Ts("Ts","tti signal",0.02,0.00000005,0.2);
  RooRealVar bb2s("bb2s","bibi2 signal",0.02, 0.0005, 0.1);
  RooRealVar fexps("fexps","f_exp signal",0.02, 0.0, 1.0);

  RooRealVar mv("mv","emme signal", 110.,10., 6000.,"GeV/c");
  RooRealVar nv("nv","enne signal", 0.93, 0.5, 15.);
  RooRealVar n2v("n2v","enne2 signal", 0.75, 0.5, 15.); 
  RooRealVar bbv("bbv","bibi signal",0.02, 0.0005, 0.1);
  RooRealVar Tv("Tv","tti signal",0.02,0.00000005,0.2);
  RooRealVar bb2v("bb2v","bibi2 signal",0.02, 0.0005, 0.1);
  RooRealVar fexpv("fexpv","f_exp signal",0.02, 0.0, 1.0);
  if (isVBFsignal) {
    mv.setVal(653.8);   mv.setConstant(kTRUE);
    nv.setVal(1.1705);   nv.setConstant(kTRUE);
    n2v.setVal(4.086);   n2v.setConstant(kTRUE);
    bbv.setVal(0.0294); // bbv.setConstant(kTRUE);
    Tv.setVal(0.0064);   // Tv.setConstant(kTRUE);
    bb2v.setVal(0.020);   // bb2v.setConstant(kTRUE);
    fexpv.setVal(0.1);   fexpv.setConstant(kTRUE);
  }

  ms.setVal(1803.8);   ms.setConstant(kTRUE);
  ns.setVal(0.733);   ns.setConstant(kTRUE);
  n2s.setVal(0.95);   n2s.setConstant(kTRUE);
  bbs.setVal(0.0323); // bbs.setConstant(kTRUE);
  Ts.setVal(0.064);   // Ts.setConstant(kTRUE);
  bb2s.setVal(0.020);   // bb2s.setConstant(kTRUE);
  fexps.setVal(0.1);   // fexps.setConstant(kTRUE);
 
  
  RooTsallis3* rt = new RooTsallis3("rt","rt",*pts,mv,nv,n2v,bbv,bb2v,Tv,fexpv);
  ws->import(*rt);
  RooTsallis3* rt2 = new RooTsallis3("rt2","rt2",*pts,ms,ns,n2s,bbs,bb2s,Ts,fexps);
  ws->import(*rt2);
  RooTsallis3* rt3 = new RooTsallis3("rt3","rt3",*pt,m,n,n2,bb,bb2,T,fexp);
  ws->import(*rt3);
  /* ws->factory("Gaussian::gau(pt,mean[3.,2.,15.],sigma[2.,0.1,10.])");
  ws->factory("Chebychev::che(pt,{a0[0.5,0.,1.],a1[0.5,0.,1.]})");
  // ws->factory("SUM::all(coeffPol[0.1,0.,1.]*che,coeffGau[0.1,0.,1.]*gau,rt)");
  ws->factory("SUM::all(coeffGau[0.1,0.,1.]*gau,rt2)");*/

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
  RooFitResult* sigVBFfit = new RooFitResult();
  if (isVBFsignal) {
    sigVBFfit = ws->pdf("rt")->fitTo(*sigVBF,Minos(0),Save(1),SumW2Error(kTRUE),NumCPU(1));
  }

  RooPlot *framesig = pts->frame();
  sprintf(reducestr,"pts > %f && pts < %f",pts->getMin(),pts->getMax());

  if (writeWeightHisto && !isVBFsignal) gg->plotOn(framesig,DataError(RooAbsData::SumW2),Cut(reducestr),LineColor(kRed),MarkerColor(kRed));
  if (!isVBFsignal) {
    sig->plotOn(framesig,DataError(RooAbsData::SumW2),Cut(reducestr));
    // ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue),Normalization(sig->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue));
  } else {
    sigVBF->plotOn(framesig,DataError(RooAbsData::SumW2),Cut(reducestr));
    // ws->pdf("rt2")->plotOn(framesig,LineColor(kBlue),Normalization(sig->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("rt")->plotOn(framesig,LineColor(kBlue));
  }
  RooHist *hpullsig = framesig->pullHist();
  RooHist *hresidsig = framesig->residHist();

  if (isVBFsignal) sprintf(fileToSave,"text/paramsVbf_%dGeV_%dTeV_all.txt",int(mZZcenter),LHCsqrts); 
  else sprintf(fileToSave,"text/paramsSig_%dGeV_%dTeV_all.txt",int(mZZcenter),LHCsqrts);

  ofstream os2(fileToSave);
  (RooArgSet(sigfit->floatParsFinal())).writeToStream(os2,false);
  if (isVBFsignal) {
    (RooArgSet(sigVBFfit->floatParsFinal())).writeToStream(os2,false);
  }
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
    float thisweight;
    double thiserrup, thiserrdown;
    if (!isVBFsignal) {
      thisweight = sig->weight(RooArgSet(*pts));
      sig->weightError(thiserrup, thiserrdown, RooAbsData::SumW2);
    } else {
      thisweight = sigVBF->weight(RooArgSet(*pts));
      sigVBF->weightError(thiserrup, thiserrdown, RooAbsData::SumW2);
    }
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
  
  RooTsallis3* rt4 = new RooTsallis3("rt4","rt4",*ws->var("pts"),*ws->var("m"),*ws->var("n"),*ws->var("n2"),*ws->var("bb"), 
                          *ws->var("bb2"),*ws->var("T"),*ws->var("fexp"));
  ws->import(*rt4);


  if (run2D) Run2D(ws,int(mZZcenter));

  /* if (writeWorkspace) {
    sprintf(fileToSave,"workspacePt_%dGeV_%dTeV_all.root",int(mZZcenter),LHCsqrts);
    TFile fout(fileToSave,"RECREATE");
    ws->Write();
    } */
  

  if (writeWeightHisto) {

    sprintf(fileToOpen,"weights/gg125_weightsUE.root");
    if (mZZcenter > 190.) sprintf(fileToOpen,"weights/gg200_weightsUE.root");
    if (mZZcenter > 390.) sprintf(fileToOpen,"weights/gg400_weightsUE.root");
    TFile fue(fileToOpen);
    TH1F* wei = (TH1F*)fue.Get("wei");

    sprintf(fileToSave,"weights/weightHisto_%dGeV_%dTeV_%s.root",int(mZZcenter),LHCsqrts,typeSyst.c_str());
    TFile fout(fileToSave,"RECREATE");
    ws->var("pts")->setMin(0.);
    ws->var("pts")->setMax(500.);

    TH1F* ggHRes = (TH1F*)ggH->Clone();
    ggH->SetName("ggH");   ggH->SetTitle("POWHEG gg pT histo");
    if (isVBFsignal) {
      sigVBFH->SetName("sigVBFH");   sigVBFH->SetTitle("POWHEG VBF pT histo");
    }
    ggHRes->SetName("ggHRes");   ggHRes->SetTitle("HRes gg pT fit histo");
    ggH->Sumw2();
    ggHRes->Sumw2();
    ggH->Scale(1./ggH->Integral());

    TH1F* wH = (TH1F*)ggH->Clone();
    wH->SetName("wH");   wH->SetTitle("HRes to POWHEG pT weight histo");  
    for (Int_t i=1; i<=ggH->GetNbinsX(); i++) {
      ws->var("pts")->setVal(ggH->GetXaxis()->GetBinCenter(i));
      ggHRes->SetBinContent(i,ws->pdf("rt2")->getVal());
      ggHRes->SetBinError(i,0.00000001);
      // cout << i << " " << ws->var("pts")->getVal() << " " << ggH->GetBinContent(i) << " " << ggHRes->GetBinContent(i) << endl;
    }
    ggHRes->Scale(1./ggHRes->Integral());

    // MODIFY SPECTRUM
    // 1) Smooth low pt?
    TH1F* ggHResMod = (TH1F*)ggHRes->Clone();
    TH1F* ggHMod = (TH1F*)ggH->Clone();
    float totalFirstBins = 0.;
    float totalFirstBinsRes = 0.;
    int nBinsToSmooth = 0;  // no Smoothing
    for (Int_t i=1; i<=nBinsToSmooth; i++) {
      totalFirstBins += ggH->GetBinContent(i);
      totalFirstBinsRes += ggHRes->GetBinContent(i);
    }
    for (Int_t i=1; i<=nBinsToSmooth; i++) {
      ggHMod->SetBinContent(i,totalFirstBins/nBinsToSmooth);
      ggHResMod->SetBinContent(i,totalFirstBinsRes/nBinsToSmooth);
    }
    
    // 2) Reweight for PS+UE/PS only (suggested in YR2)
    ggHResMod->Multiply(ggHResMod,wei);

    wH->Divide(ggHResMod,ggHMod);
    ggHMod->Write();
    ggHResMod->Write();
    if (isVBFsignal) sigVBFH->Write();
    wH->Write();
    fout.Close();
  }

}



