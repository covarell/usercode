// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>

void studyPtSyst(int whichtype = 1, bool also7TeV = true) {

  // -5 - PDF : VBF
  // -4 - scales : VBF
  // -3 - effect of finite top mass
  // -2 - mu_Q in gg signal shape
  // -1 - fraction of VBF (NOT for VBF fraction fitting)
  // 1 - Z+X in background
  // 2 - Cross-check Z+X vs. Z+jets
  // 3 - low pT shape in background: nonblinded
  // 4 - low pT shape in background: single Z
  // 5 - adding/not adding ggZZ
  // 6 - PDF : ZZ
  // 7 - scales : ZZ

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
 
  // READ TREES
  TChain* ggTree = new TChain("angles");
  if (also7TeV) {
    ggTree->Add("../datafiles/4e/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/4mu/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/2mu2e/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
  } else {
    ggTree->Add("../datafiles/4e/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/4mu/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/2mu2e/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root"); 
  }

  TChain* VBFTree = new TChain("angles");
  if (also7TeV) {
    VBFTree->Add("../datafiles/4e/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/4mu/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/2mu2e/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
  } else {
    VBFTree->Add("../datafiles/4e/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/4mu/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/2mu2e/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root"); 
  }

  TChain* zzTree = new TChain("angles");
  if (also7TeV) {
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4e_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo2e2tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4mu_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo2mu2tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2mu_105-1000_7TeV_withDiscriminantsPtY.root"); 
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2mu2tau_105-1000_7TeV_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo4tau_105-1000_7TeV_withDiscriminantsPtY.root");
  } else {
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4e_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo2e2tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4mu_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo2mu2tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2mu_105-1000_withDiscriminantsPtY.root"); 
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2mu2tau_105-1000_withDiscriminantsPtY.root");
    zzTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo4tau_105-1000_withDiscriminantsPtY.root");
  }

  TChain* ggzzTree = new TChain("angles");
  if (also7TeV) {
    ggzzTree->Add("../datafiles/4e/HZZ4lTree_ggZZ4l_105-1000_7TeV_withDiscriminantsPtY.root");
    ggzzTree->Add("../datafiles/4mu/HZZ4lTree_ggZZ4l_105-1000_7TeV_withDiscriminantsPtY.root");
    ggzzTree->Add("../datafiles/2mu2e/HZZ4lTree_ggZZ2l2l_105-1000_7TeV_withDiscriminantsPtY.root");
  } else {
    ggzzTree->Add("../datafiles/4e/HZZ4lTree_ggZZ4l_105-1000_withDiscriminantsPtY.root");
    ggzzTree->Add("../datafiles/4mu/HZZ4lTree_ggZZ4l_105-1000_withDiscriminantsPtY.root");
    ggzzTree->Add("../datafiles/2mu2e/HZZ4lTree_ggZZ2l2l_105-1000_withDiscriminantsPtY.root");
  }
   
  TChain* dataTree = new TChain("angles");
  if (also7TeV) {
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleEle_105-1000_7TeV_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleMu_105-1000_7TeV_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleOr_105-1000_7TeV_withDiscriminantsPtY.root");
  } else {
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleEle_105-1000_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleMu_105-1000_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleOr_105-1000_withDiscriminantsPtY.root"); 
  }

  TChain* crTree = new TChain("angles");
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEEEssTree_105-1000_withDiscriminantsPtY.root");
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEMMssTree_105-1000_withDiscriminantsPtY.root");
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEEEssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEMMssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEEEssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEMMssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");	
  crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMMMssTree_105-1000_withDiscriminantsPtY.root"); 
  if (also7TeV) {
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMEEssTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMMMssTree_105-1000_7TeV_withDiscriminantsPtY.root");
  }

  TChain* crosTree = new TChain("angles");
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEEEosTree_105-1000_withDiscriminantsPtY.root");
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEMMosTree_105-1000_withDiscriminantsPtY.root");
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMEEosTree_105-1000_withDiscriminantsPtY.root");
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMMMosTree_105-1000_withDiscriminantsPtY.root");
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEEEosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEMMosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMEEosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMMMosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEEEosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEMMosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMEEosTree_105-1000_withDiscriminantsPtY.root");	
  crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMMMosTree_105-1000_withDiscriminantsPtY.root");	
  if (also7TeV) {
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CREEMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleEle_CRMMMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CREEMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleMu_CRMMMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CREEMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMEEosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
    crosTree->Add("../datafiles/CR/HZZ4lTree_DoubleOr_CRMMMMosTree_105-1000_7TeV_withDiscriminantsPtY.root");	
  }

  TChain* crzjTree = new TChain("angles");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEEEssTree_105-1000_withDiscriminantsPtY.root");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEMMssTree_105-1000_withDiscriminantsPtY.root");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEMMssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEMMssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEMMssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMEEssTree_105-1000_withDiscriminantsPtY.root");	
  crzjTree->Add("../datafiles/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMMMssTree_105-1000_withDiscriminantsPtY.root");

  float mgg, mVBF, mzz, mggzz, mdata, mcr, mcros, mzj;
  float wgg, wVBF, wzz, wggzz, wzj;
  float nlogg, nloVBF, nlozz, nloggzz, nlodata, nlocr, nlozj; 
  float ptgg, ptVBF, ptzz, ptggzz, ptdata, ptcr, ptzj;

  ggTree->SetBranchAddress("zzmass",&mgg);
  ggTree->SetBranchAddress("MC_weight",&wgg);
  ggTree->SetBranchAddress("ZZPt",&ptgg);
  ggTree->SetBranchAddress("melaLDWithPtY",&nlogg);

  VBFTree->SetBranchAddress("zzmass",&mVBF);
  VBFTree->SetBranchAddress("MC_weight",&wVBF);
  VBFTree->SetBranchAddress("ZZPt",&ptVBF);
  VBFTree->SetBranchAddress("melaLDWithPtY",&nloVBF);

  zzTree->SetBranchAddress("zzmass",&mzz);
  zzTree->SetBranchAddress("MC_weight",&wzz);
  zzTree->SetBranchAddress("ZZPt",&ptzz);
  zzTree->SetBranchAddress("melaLDWithPtY",&nlozz);

  ggzzTree->SetBranchAddress("zzmass",&mggzz);
  ggzzTree->SetBranchAddress("MC_weight",&wggzz);
  ggzzTree->SetBranchAddress("ZZPt",&ptggzz);
  ggzzTree->SetBranchAddress("melaLDWithPtY",&nloggzz);

  dataTree->SetBranchAddress("zzmass",&mdata);
  dataTree->SetBranchAddress("ZZPt",&ptdata);
  dataTree->SetBranchAddress("melaLDWithPtY",&nlodata);

  crTree->SetBranchAddress("zzmass",&mcr);
  crTree->SetBranchAddress("ZZPt",&ptcr);
  crTree->SetBranchAddress("melaLDWithPtY",&nlocr);

  crosTree->SetBranchAddress("zzmass",&mcros);

  crzjTree->SetBranchAddress("zzmass",&mzj);
  crzjTree->SetBranchAddress("MC_weight",&wzj);
  crzjTree->SetBranchAddress("ZZPt",&ptzj);
  crzjTree->SetBranchAddress("melaLDWithPtY",&nlozj);

  char nameSyst[200]; 
  char nameFile[200];

  sprintf(nameFile,"../../weightHisto_125GeV_8TeV_all.root");
  if (also7TeV) sprintf(nameFile,"../../weightHisto_125GeV_7TeV_all.root");
  TFile weightsSig(nameFile);
  TH1F* wHalf = (TH1F*)weightsSig.Get("wH");
  TFile weightsSig1("../../weightHisto_125GeV_8TeV_One.root");
  TH1F* wOne = (TH1F*)weightsSig1.Get("wH");
  TFile weightsSig2("../../weightHisto_125GeV_8TeV_Two.root");
  TH1F* wTwo = (TH1F*)weightsSig2.Get("wH");

  TFile pwhgNoMass("ggH125_infiniteMT.root");
  TFile pwhgMass("ggH125_finiteMT.root");
    
  TH1F* nmH = (TH1F*)((TH2F*)pwhgNoMass.Get("Pt_sig"))->ProjectionY("nmH");
  TH1F* mH = (TH1F*)((TH2F*)pwhgMass.Get("Pt_sig"))->ProjectionY("mH");
 
  nmH->Rebin(4);   nmH->Sumw2();
  mH->Rebin(4);    mH->Sumw2();
  nmH->GetXaxis()->SetRange(1,80);
  mH->GetXaxis()->SetRange(1,80);
  nmH->Scale(1./nmH->Integral());
  mH->Scale(1./mH->Integral());
  TH1F* ratio = (TH1F*)mH->Clone();  
  ratio->Divide(mH,nmH);

  TCanvas can("can","The canvas",5.,5.,600.,700.); 
  can.Divide(1,2);

  int nbins2 = 160;
  float theMax = 400.;

  // Standard pT binning
  // TH1F* pth = new TH1F("pth","pt",50,0.,200.);
  TH1F* pth1 = new TH1F("pth1","pt",nbins2,0.,theMax);
  // TH1F* pth2 = new TH1F("pth2","pt",50,0.,200.);
  //Extended ranges
  TH1F* pth = new TH1F("pth","pt",nbins2,0.,theMax);
  TH1F* pth2 = new TH1F("pth2","pt",nbins2,0.,theMax);
  TH1F* ptvbf = new TH1F("ptvbf","pt",nbins2,0.,theMax);
  pth->Sumw2();
  pth1->Sumw2();
  pth2->Sumw2();
  ptvbf->Sumw2();
 
  //Alternative pT binnings
  const int nbins = 10;
  float binlimits[nbins+1] = {0.,10.,20.,30.,40.,50.,70.,90.,110.,150.,190.}; 
  TH1F* ptalt = new TH1F("ptalt","pt",nbins,binlimits);
  TH1F* ptalt1 = new TH1F("ptalt1","pt",nbins,binlimits);
  ptalt->Sumw2();
  ptalt1->Sumw2();
  
  TH1F* nloh = new TH1F("nloh","nloMELA",52,-0.02,1.02);
  TH1F* nloh1 = new TH1F("nloh1","nloMELA",52,-0.02,1.02);
  TH1F* nloh2 = new TH1F("nloh2","nloMELA",52,-0.02,1.02);
  nloh->Sumw2();
  nloh1->Sumw2();
  nloh2->Sumw2();

  if (whichtype == -5) {

    TFile* pdfsfile[3];
    pdfsfile[0] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_standard.root");
    pdfsfile[1] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_MSTW.root");
    pdfsfile[2] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_NNPDF.root");

    TH1F* pdfh[3];
    TH1F* ratiopdf[2];
    for (int i = 0; i < 3; i++) {
      sprintf(nameFile,"pdfh%d",i);
      pdfh[i] = (TH1F*)((TH2F*)pdfsfile[i]->Get("Pt_sig"))->ProjectionY(nameFile);
      pdfh[i]->Rebin(4);   pdfh[i]->Sumw2();
      pdfh[i]->GetXaxis()->SetRange(1,80);
      pdfh[i]->Scale(1./pdfh[i]->Integral());
      if (i > 0) {
	ratiopdf[i-1] = (TH1F*)pdfh[i]->Clone();
	ratiopdf[i-1]->Divide(pdfh[i],pdfh[0]);
      } 
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptVBF < 500.) {
	int theBin = pdfh[0]->FindBin(ptVBF);
 	pth->Fill(ptVBF,wVBF);
	pth1->Fill(ptVBF,wVBF*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(ptVBF,wVBF*ratiopdf[0]->GetBinContent(theBin));
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*ratiopdf[1]->GetBinContent(theBin));
	nloh2->Fill(nloVBF,wVBF*ratiopdf[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVBF,wVBF*thisDiff); 
      }
    } 

    sprintf(nameFile,"PT_vbf_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_vbf_SEL_7TeV.root");
    TFile f1(nameFile,"UPDATE");
    pth->SetName("ptH_Default");
    pth->Write();
    ptvbf->SetName("ptH_PDF");
    ptvbf->Write();
    f1.Close();

    sprintf(nameSyst,"PDF-VBF");

  }

  if (whichtype == -4) {

    TFile* scalesfile[7];
    scalesfile[0] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_standard.root");
    scalesfile[1] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_mufdouble.root");
    scalesfile[2] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_mufhalf.root");
    scalesfile[3] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_murdouble.root");
    scalesfile[4] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_murhalf.root");
    scalesfile[5] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_mufdouble_murdouble.root");
    scalesfile[6] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/VBF_mufhalf_murhalf.root");

    TH1F* scaleh[7];
    TH1F* ratioscale[6];
    for (int i = 0; i < 7; i++) {
      sprintf(nameFile,"scaleh%d",i);
      scaleh[i] = (TH1F*)((TH2F*)scalesfile[i]->Get("Pt_sig"))->ProjectionY(nameFile);
      scaleh[i]->Rebin(4);   scaleh[i]->Sumw2();
      scaleh[i]->GetXaxis()->SetRange(1,80);
      scaleh[i]->Scale(1./scaleh[i]->Integral());
      if (i > 0) {
	ratioscale[i-1] = (TH1F*)scaleh[i]->Clone();
	ratioscale[i-1]->Divide(scaleh[i],scaleh[0]);
      } 
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptVBF < 500.) {
	int theBin = scaleh[0]->FindBin(ptVBF);
 	pth->Fill(ptVBF,wVBF);
	pth1->Fill(ptVBF,wVBF*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(ptVBF,wVBF*ratioscale[0]->GetBinContent(theBin));
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*ratioscale[1]->GetBinContent(theBin));
	nloh2->Fill(nloVBF,wVBF*ratioscale[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVBF,wVBF*thisDiff); 
      }
    } 

    sprintf(nameFile,"PT_vbf_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_vbf_SEL_7TeV.root");
    TFile f1(nameFile,"UPDATE");
    ptvbf->SetName("ptH_scale");
    ptvbf->Write();
    f1.Close();

    sprintf(nameSyst,"scale-VBF");

  }  

  if (whichtype == -3) {

    // Grazzini's matching <-- NO! (Vicini)
    /* for (Int_t iBin = 1; iBin <= nmH->GetNbinsX() ; ++iBin) {
      float thept = nmH->GetBinCenter(iBin);
      float matchw = 1.;
      // float matchw = 0.;   // should be zero, but allow division
      // if (thept > 149.8) matchw = 1.;
      //if (thept < 149.8 && thept > 89.8) 
      // matchw = 0.5*(1 - cos(3.1416*(thept-89.8)/60.));
      // cout << matchw << endl;
      nmH->SetBinContent(iBin,nmH->GetBinContent(iBin)*matchw);
      mH->SetBinContent(iBin,mH->GetBinContent(iBin)*matchw);
      } */

    // nmH->Scale(1./nmH->Integral());
    // mH->Scale(1./mH->Integral());

    mH->SetMarkerColor(1);    mH->SetLineColor(1);    mH->SetLineWidth(2);
    nmH->SetMarkerColor(2);   nmH->SetLineColor(2);   nmH->SetLineWidth(2);

    can.cd(1);
    gPad->SetBottomMargin(0.0);
    mH->GetXaxis()->SetLabelColor(kWhite);
    mH->GetYaxis()->SetTitle("Normalized entries");
    mH->Draw();
    nmH->Draw("SAME");

    TH1F* diffpt = (TH1F*)mH->Clone();
    TH1F* pullpt = (TH1F*)mH->Clone();
    diffpt->Add(nmH,mH,1,-1);
    pullpt->Divide(diffpt,nmH);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-1.);
    pullpt->SetMaximum(1.);
    pullpt->GetXaxis()->SetLabelColor(kBlack);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");

    can.SaveAs("pt_systTopMass_corr.pdf");
    sprintf(nameSyst,"TopMass");

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
        float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptgg,wgg*wHalf->GetBinContent(theBin));
	pth2->Fill(ptgg,wgg*newW*wHalf->GetBinContent(theBin));
	ptvbf->Fill(ptgg,wgg*newW*wHalf->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*wHalf->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
      }
    } 

    sprintf(nameFile,"PT_gg_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_gg_SEL_7TeV.root");
    TFile f(nameFile,"UPDATE");
    pth2->SetName("ptH_TopMass");
    pth2->Write();
    f.Close();

    /* for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth->Fill(ptVBF,wVBF);
	pth2->Fill(ptVBF,wVBF);
        nloh->Fill(nloVBF,wVBF);
	nloh2->Fill(nloVBF,wVBF);
      }
      }*/

    // return;

  }  

  if (whichtype == -2) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptgg,wgg*newW*wHalf->GetBinContent(theBin));
	pth1->Fill(ptgg,wgg*newW*wOne->GetBinContent(theBin));
	pth2->Fill(ptgg,wgg*newW*wTwo->GetBinContent(theBin));
	ptvbf->Fill(ptgg,wgg*newW*wTwo->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*newW*wOne->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*wTwo->GetBinContent(theBin));
      }
    } 

    sprintf(nameFile,"PT_gg_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_gg_SEL_7TeV.root");
    TFile f(nameFile,"UPDATE");
    pth->SetName("ptH_Default");
    pth->Write();
    pth2->SetName("ptH_Resummation");
    pth2->Write();
    f.Close();

    /* for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth->Fill(ptVBF,wVBF);
        pth1->Fill(ptVBF,wVBF);
	pth2->Fill(ptVBF,wVBF);
	ptvbf->Fill(ptVBF,wVBF);
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF);
	nloh2->Fill(nloVBF,wVBF);
      }
      }*/

    sprintf(nameSyst,"Resummation");
 
  }

  if (whichtype == -1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptgg,wgg*newW*wHalf->GetBinContent(theBin));
	pth1->Fill(ptgg,wgg*newW*0.84*wHalf->GetBinContent(theBin));  // gg: +/-16%
	pth2->Fill(ptgg,wgg*newW*1.16*wHalf->GetBinContent(theBin));
	ptvbf->Fill(ptgg,wgg*newW*1.16*wHalf->GetBinContent(theBin));  
        nloh->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*newW*0.84*wHalf->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*1.16*wHalf->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth->Fill(ptVBF,wVBF);
        pth1->Fill(ptVBF,wVBF*1.02);             // VBF: +/-2%
	pth2->Fill(ptVBF,wVBF*0.98);
	ptvbf->Fill(ptVBF,wVBF*0.98);
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*1.02);
	nloh2->Fill(nloVBF,wVBF*0.98);
      }
    }

    sprintf(nameSyst,"VBFfraction");
 
  }


  if (whichtype == 1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth1->Fill(ptgg,wgg*newW*wHalf->GetBinContent(theBin));
        nloh1->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth1->Fill(ptVBF,wVBF);
        nloh1->Fill(nloVBF,wVBF);
      }
    }

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	pth2->Fill(ptzz,wzz);
	ptvbf->Fill(ptzz,wzz);
        nloh2->Fill(nlozz,wzz);
      }
    } 

    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f2(nameFile,"UPDATE");
    pth2->SetName("ptH_Default");
    pth2->Write();
    f2.Close();

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < 140. && mcr > 105. && ptcr < 500.) {
	pth->Fill(ptcr);
        nloh->Fill(nlocr);
      }
    } 

    sprintf(nameFile,"PT_zx_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zx_SEL_7TeV.root");
    TFile f3(nameFile,"UPDATE");
    pth->SetName("ptH_Default");
    pth->Write();
    f3.Close();

    sprintf(nameSyst,"ZplusX");
 
  }

  if (whichtype == 2) {

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (ptcr < 500.) {          // do it in the whole mass range
	pth->Fill(ptcr);
        nloh->Fill(nlocr);
      }
    } 

    for (Int_t iEvt = 0; iEvt < crzjTree->GetEntries() ; ++iEvt) {
      crzjTree->GetEntry(iEvt);
      if (ptzj < 500.) {
	pth2->Fill(ptzj,wzj);
	ptvbf->Fill(ptzj,wzj);
        nloh2->Fill(nlozj,wzj);	
      }
    }

    sprintf(nameSyst,"ZjetsXcheck");
 
  }

  if (whichtype == 3) {

    sprintf(nameSyst,"histograms8TeV.root");
    if (also7TeV) sprintf(nameSyst,"histograms7TeV.root");
    TFile fileMike(nameSyst);

    sprintf(nameSyst,"tmpTH1");
    if (also7TeV) sprintf(nameSyst,"DATA");
    TH1F* dataUnb = (TH1F*)fileMike.Get(nameSyst);
    dataUnb->GetXaxis()->SetRange(1,20);
    TH1F* zzUnb = (TH1F*)fileMike.Get("ZZ");
    zzUnb->Scale(13.6/zzUnb->GetBinContent(4));
    zzUnb->GetXaxis()->SetRange(1,20);
    TH1F* zxUnb = (TH1F*)fileMike.Get("FAKES");
    zxUnb->Scale(1.7/zxUnb->GetBinContent(4));
    zxUnb->GetXaxis()->SetRange(1,20);

    // Subtract fakes 
    dataUnb->Add(dataUnb,zxUnb,1,-1);
    dataUnb->SetMarkerColor(1);    dataUnb->SetLineColor(1);    dataUnb->SetLineWidth(2);
    zzUnb->SetMarkerColor(2);   zzUnb->SetLineColor(2);   zzUnb->SetLineWidth(2);
    can.cd(1);
    gPad->SetBottomMargin(0.0);
    dataUnb->GetXaxis()->SetLabelColor(kWhite);
    dataUnb->GetYaxis()->SetTitle("Entries");
    dataUnb->Draw();
    zzUnb->Draw("SAME");

    TH1F* diffpt = (TH1F*)zxUnb->Clone();
    TH1F* pullpt = (TH1F*)zxUnb->Clone();
    diffpt->Add(dataUnb,zzUnb,1,-1);
    pullpt->Divide(diffpt,dataUnb);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-2.);
    pullpt->SetMaximum(2.);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
    TF1 *mypol1 = new TF1("mypol1","pol1");
    pullpt->Fit("mypol1","","",0.,110.);

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	pth->Fill(ptzz,wzz*(1. + mypol1->Eval(ptzz)));
      }
    } 
    
    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f2(nameFile,"UPDATE");
    pth->SetName("ptH_UnbRegion");
    pth->Write();
    f2.Close();

    can.SaveAs("pt_systUnbRegion.pdf");
    return;

  }

  if (whichtype == 4) {
    
    float dataSz[nbins] = {465.0,251.5,116.0,59.8,33.8,18.1,7.79,4.75,1.93,0.60};
    float dataErr[nbins] = {7.5,5.4,4.0,2.7,1.8,0.9,0.54,0.42,0.17,0.10};
    float powSz[nbins] = {454.3,256.7,118.0,58.0,32.0,18.2,7.6,3.4,1.2,0.59};
    float powErr[nbins] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    for (Int_t iBin = 0; iBin < nbins; ++iBin) {
      ptalt->SetBinContent(iBin+1,dataSz[iBin]);
      ptalt->SetBinError(iBin+1,dataErr[iBin]);
      ptalt1->SetBinContent(iBin+1,powSz[iBin]);
      ptalt1->SetBinError(iBin+1,powErr[iBin]);
    }

    ptalt->SetMarkerColor(1);    ptalt->SetLineColor(1);    ptalt->SetLineWidth(2);
    ptalt1->SetMarkerColor(2);   ptalt1->SetLineColor(2);   ptalt1->SetLineWidth(2);
    can.cd(1);
    gPad->SetBottomMargin(0.0);
    ptalt->GetXaxis()->SetLabelColor(kWhite);
    ptalt->GetYaxis()->SetTitle("Entries");
    ptalt->Draw();
    ptalt1->Draw("SAME");

    TH1F* diffpt = (TH1F*)ptalt->Clone();
    TH1F* pullpt = (TH1F*)ptalt->Clone();
    diffpt->Add(ptalt,ptalt1,1,-1);
    pullpt->Divide(diffpt,ptalt);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-1.);
    pullpt->SetMaximum(1.);
    pullpt->GetXaxis()->SetLabelColor(kBlack);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
    TF1 *mypol1 = new TF1("mypol1","pol1");
    pullpt->Fit("mypol1","","");

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	pth->Fill(ptzz,wzz*(1. + mypol1->Eval(ptzz)));
      }
    } 
    
    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f2(nameFile,"UPDATE");
    pth->SetName("ptH_SingleZ");
    pth->Write();
    f2.Close();

    can.SaveAs("pt_systSingleZ.pdf");
    return;
  }

  if (whichtype == 5) {

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	pth->Fill(ptzz,wzz);
        nloh->Fill(nlozz,wzz);
	pth2->Fill(ptzz,wzz);
	ptvbf->Fill(ptzz,wzz);
        nloh2->Fill(nlozz,wzz);
      }
    } 

    for (Int_t iEvt = 0; iEvt < ggzzTree->GetEntries() ; ++iEvt) {
      ggzzTree->GetEntry(iEvt);
      if (mggzz < 140. && mggzz > 105. && ptggzz < 500.) {
	pth1->Fill(ptggzz,wggzz);
        nloh1->Fill(nloggzz,wggzz);
	pth->Fill(ptggzz,wggzz);
        nloh->Fill(nloggzz,wggzz);
      }
    } 

    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f4(nameFile,"UPDATE");
    pth->SetName("ptH_ggZZ");
    pth->Write();
    f4.Close();

    sprintf(nameSyst,"ggZZ");
  }

  if (whichtype == 6) {

    TFile* pdfsfile[3];
    pdfsfile[0] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_standard.root");
    pdfsfile[1] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_MSTW.root");
    pdfsfile[2] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_NNPDF.root");

    TH1F* pdfh[3];
    TH1F* ratiopdf[2];
    for (int i = 0; i < 3; i++) {
      sprintf(nameFile,"pdfh%d",i);
      pdfh[i] = (TH1F*)((TH2F*)pdfsfile[i]->Get("Pt_bkg"))->ProjectionY(nameFile);
      pdfh[i]->Rebin(4);   pdfh[i]->Sumw2();
      pdfh[i]->GetXaxis()->SetRange(1,80);
      pdfh[i]->Scale(1./pdfh[i]->Integral());
      if (i > 0) {
	ratiopdf[i-1] = (TH1F*)pdfh[i]->Clone();
	ratiopdf[i-1]->Divide(pdfh[i],pdfh[0]);
      } 
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	int theBin = pdfh[0]->FindBin(ptzz);
 	pth->Fill(ptzz,wzz);
	pth1->Fill(ptzz,wzz*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(ptzz,wzz*ratiopdf[0]->GetBinContent(theBin));
        nloh->Fill(nlozz,wzz);
	nloh1->Fill(nlozz,wzz*ratiopdf[1]->GetBinContent(theBin));
	nloh2->Fill(nlozz,wzz*ratiopdf[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptzz,wzz*thisDiff); 
      }
    } 

    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f1(nameFile,"UPDATE");
    ptvbf->SetName("ptH_PDF");
    ptvbf->Write();
    f1.Close();

    sprintf(nameSyst,"PDF-ZZ");

  }

  if (whichtype == 7) {

    TFile* scalesfile[7];
    scalesfile[0] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_standard.root");
    scalesfile[1] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_mufdouble.root");
    scalesfile[2] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_mufhalf.root");
    scalesfile[3] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_murdouble.root");
    scalesfile[4] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_murhalf.root");
    scalesfile[5] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_mufdouble_murdouble.root");
    scalesfile[6] = TFile::Open("../../../../src/GeneratorInterface/ExternalDecays/test/ZZ_mufhalf_murhalf.root");

    TH1F* scaleh[7];
    TH1F* ratioscale[6];
    for (int i = 0; i < 7; i++) {
      sprintf(nameFile,"scaleh%d",i);
      scaleh[i] = (TH1F*)((TH2F*)scalesfile[i]->Get("Pt_bkg"))->ProjectionY(nameFile);
      scaleh[i]->Rebin(4);   scaleh[i]->Sumw2();
      scaleh[i]->GetXaxis()->SetRange(1,80);
      scaleh[i]->Scale(1./scaleh[i]->Integral());
      if (i > 0) {
	ratioscale[i-1] = (TH1F*)scaleh[i]->Clone();
	ratioscale[i-1]->Divide(scaleh[i],scaleh[0]);
      } 
    }

    // Draw and take the maximum difference per bin
    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 140. && mzz > 105. && ptzz < 500.) {
	int theBin = scaleh[0]->FindBin(ptzz);
 	pth->Fill(ptzz,wzz);
	pth1->Fill(ptzz,wzz*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(ptzz,wzz*ratioscale[0]->GetBinContent(theBin));
        nloh->Fill(nlozz,wzz);
	nloh1->Fill(nlozz,wzz*ratioscale[1]->GetBinContent(theBin));
	nloh2->Fill(nlozz,wzz*ratioscale[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptzz,wzz*thisDiff); 
      }
    } 

    sprintf(nameFile,"PT_zz_SEL_8TeV.root");
    if (also7TeV) sprintf(nameFile,"PT_zz_SEL_7TeV.root");
    TFile f1(nameFile,"UPDATE");
    ptvbf->SetName("ptH_scale");
    ptvbf->Write();
    f1.Close();

    sprintf(nameSyst,"scale-ZZ");

  }

  pth->Scale(1./pth->Integral());
  pth1->Scale(1./pth1->Integral());
  pth2->Scale(1./pth2->Integral());
  ptvbf->Scale(1./ptvbf->Integral());
  nloh->Scale(1./nloh->Integral());
  nloh1->Scale(1./nloh1->Integral());
  nloh2->Scale(1./nloh2->Integral());
  
  pth->SetMarkerColor(1);    pth->SetLineColor(1);    pth->SetLineWidth(2);
  pth1->SetMarkerColor(2);   pth1->SetLineColor(2);   pth1->SetLineWidth(2);
  pth2->SetMarkerColor(4);   pth2->SetLineColor(4);   pth2->SetLineWidth(2);
  			    			   
  nloh->SetMarkerColor(1);   nloh->SetLineColor(1);   nloh->SetLineWidth(2);
  nloh1->SetMarkerColor(2);  nloh1->SetLineColor(2);  nloh1->SetLineWidth(2);
  nloh2->SetMarkerColor(4);  nloh2->SetLineColor(4);  nloh2->SetLineWidth(2);
  
  can.cd(1);
  gPad->SetBottomMargin(0.0);
  if (whichtype == -3) gPad->SetLogy();
  pth2->GetXaxis()->SetLabelColor(kWhite);
  pth2->GetYaxis()->SetTitle("Normalized entries");
  pth2->Draw();
  if (whichtype != -3) pth1->Draw("SAME");
  pth->Draw("SAME");
  
  TH1F* diffpt = (TH1F*)pth->Clone();
  TH1F* pullpt = (TH1F*)pth->Clone();
  diffpt->Add(ptvbf,pth,1,-1);
  pullpt->Divide(diffpt,pth);
  can.cd(2);
  gPad->SetTopMargin(0.0);
  pullpt->SetMinimum(-1.);
  pullpt->SetMaximum(1.);
  pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  pullpt->GetYaxis()->SetTitle("Relative difference");
  pullpt->Draw("E");
  
  sprintf(nameFile,"pt_syst%s.pdf",nameSyst);
  can.SaveAs(nameFile);
  
  can.cd(1);
  gPad->SetBottomMargin(0.0);
  nloh1->GetXaxis()->SetLabelColor(kWhite);
  nloh1->GetYaxis()->SetTitle("Normalized entries");
  nloh2->GetXaxis()->SetLabelColor(kWhite);
  nloh2->GetYaxis()->SetTitle("Normalized entries");
  if (whichtype == 2) nloh2->Draw();
  else nloh->Draw();
  if (whichtype != -3) nloh1->Draw();
  if (whichtype == 2) nloh->Draw("SAME");
  nloh2->Draw("SAME");
  
  TH1F* diffnlo = (TH1F*)nloh->Clone();
  TH1F* pullnlo = (TH1F*)nloh->Clone();
  diffnlo->Add(nloh2,nloh,1,-1);
  pullnlo->Divide(diffnlo,nloh);
  can.cd(2);
  gPad->SetTopMargin(0.0);
  pullnlo->SetMinimum(-1.);
  pullnlo->SetMaximum(1.);
  pullnlo->GetXaxis()->SetTitle("nloMELA");
  pullnlo->GetYaxis()->SetTitle("Relative difference");
  pullnlo->Draw("E");
  
  sprintf(nameFile,"nlo_syst%s.pdf",nameSyst);
  can.SaveAs(nameFile);
  
  if (whichtype == 1) {

    TH1F* balpt = (TH1F*)pth->Clone();
    balpt->Add(pth2,pth,20.516,4.193);
    balpt->Scale(1./balpt->Integral());

    can.cd(1);
    gPad->SetBottomMargin(0.0);
    balpt->GetXaxis()->SetLabelColor(kWhite);
    balpt->GetYaxis()->SetTitle("Normalized entries");
    pth2->Draw();
    balpt->Draw("SAME");
    pth1->Draw("SAME");

    diffpt->Add(balpt,pth2,1,-1);
    pullpt->Divide(diffpt,pth2);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-2.);
    pullpt->SetMaximum(2.);
    pullpt->Draw("E");
       
    sprintf(nameFile,"pt_balanced%s.pdf",nameSyst);
    can.SaveAs(nameFile);
    
    TH1F* balnlo = (TH1F*)nloh->Clone();
    balnlo->Add(nloh2,nloh,20.516,4.193);
    balnlo->Scale(1./balnlo->Integral());

    can.cd(1);
    gPad->SetBottomMargin(0.0);
    balnlo->GetXaxis()->SetLabelColor(kWhite);
    balnlo->GetYaxis()->SetTitle("Normalized entries");
    nloh1->Draw();
    balnlo->Draw("SAME");
    nloh2->Draw("SAME");
    
    diffnlo->Add(balnlo,nloh2,1,-1);
    pullnlo->Divide(diffnlo,nloh2);
    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullnlo->SetMinimum(-1.);
    pullnlo->SetMaximum(1.);
    pullnlo->Draw("E");

    sprintf(nameFile,"nlo_balanced%s.pdf",nameSyst);
    can.SaveAs(nameFile);
  }

  return;
}

