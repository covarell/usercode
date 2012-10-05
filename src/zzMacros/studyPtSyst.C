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

  // -2 - mu_Q in gg signal shape
  // -1 - fraction of VBF (NOT for VBF fraction fitting)
  // 1 - Z+X in background
  // 2 - Cross-check Z+X vs. Z+jets
  // 3 - low pT shape in background: nonblinded
  // 4 - low pT shape in background: single Z

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
 
  // READ TREES
  TChain* ggTree = new TChain("angles");
  ggTree->Add("../datafiles/4e/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root");
  ggTree->Add("../datafiles/4mu/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root");
  ggTree->Add("../datafiles/2mu2e/HZZ4lTree_H125_105-1000_withDiscriminantsPtY.root"); 
  if (also7TeV) {
    ggTree->Add("../datafiles/4e/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/4mu/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
    ggTree->Add("../datafiles/2mu2e/HZZ4lTree_H125_105-1000_7TeV_withDiscriminantsPtY.root");
  }

  TChain* VBFTree = new TChain("angles");
  VBFTree->Add("../datafiles/4e/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root");
  VBFTree->Add("../datafiles/4mu/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root");
  VBFTree->Add("../datafiles/2mu2e/HZZ4lTree_VBFH125_105-1000_withDiscriminantsPtY.root"); 
  if (also7TeV) {
    VBFTree->Add("../datafiles/4e/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/4mu/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
    VBFTree->Add("../datafiles/2mu2e/HZZ4lTree_VBFH125_105-1000_7TeV_withDiscriminantsPtY.root");
  }

  TChain* zzTree = new TChain("angles");
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
  } 

  TChain* dataTree = new TChain("angles");
  dataTree->Add("../datafiles/data/HZZ4lTree_DoubleEle_105-1000_withDiscriminantsPtY.root");
  dataTree->Add("../datafiles/data/HZZ4lTree_DoubleMu_105-1000_withDiscriminantsPtY.root");
  dataTree->Add("../datafiles/data/HZZ4lTree_DoubleOr_105-1000_withDiscriminantsPtY.root"); 
  if (also7TeV) {
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleEle_105-1000_7TeV_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleMu_105-1000_7TeV_withDiscriminantsPtY.root");
    dataTree->Add("../datafiles/data/HZZ4lTree_DoubleOr_105-1000_7TeV_withDiscriminantsPtY.root");
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

  float mgg, mVBF, mzz, mdata, mcr, mcros, mzj;
  float wgg, wVBF, wzz, wzj;
  float nlogg, nloVBF, nlozz, nlodata, nlocr, nlozj; 
  float ptgg, ptVBF, ptzz, ptdata, ptcr, ptzj;

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


  TFile weightsSig("../../weightHisto_125GeV_8TeV_all.root");
  TH1F* wHalf = (TH1F*)weightsSig.Get("wH");
  TFile weightsSig1("../../weightHisto_125GeV_8TeV_One.root");
  TH1F* wOne = (TH1F*)weightsSig1.Get("wH");
  TFile weightsSig2("../../weightHisto_125GeV_8TeV_Two.root");
  TH1F* wTwo = (TH1F*)weightsSig2.Get("wH");

  TCanvas can("can","The canvas",5.,5.,600.,700.); 
  can.Divide(1,2);

  // Standard pT binning
  TH1F* pth = new TH1F("pth","pt",50,0.,200.);
  TH1F* pth1 = new TH1F("pth1","pt",50,0.,200.);
  TH1F* pth2 = new TH1F("pth2","pt",50,0.,200.);
  pth->Sumw2();
  pth1->Sumw2();
  pth2->Sumw2();

  //Alternative pT binning
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

  char nameSyst[200]; 
  char nameFile[200];

  if (whichtype == -2) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
	pth->Fill(ptgg,wgg*wHalf->GetBinContent(theBin));
	pth1->Fill(ptgg,wgg*wOne->GetBinContent(theBin));
	pth2->Fill(ptgg,wgg*wTwo->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*wOne->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*wTwo->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth->Fill(ptVBF,wVBF);
        pth1->Fill(ptVBF,wVBF);
	pth2->Fill(ptVBF,wVBF);
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF);
	nloh2->Fill(nloVBF,wVBF);
      }
    }

    sprintf(nameSyst,"Resummation");
 
  }

  if (whichtype == -1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < 140. && mgg > 105. && ptgg < 500.) {
	int theBin = wHalf->FindBin(ptgg);
	pth->Fill(ptgg,wgg*wHalf->GetBinContent(theBin));
	pth1->Fill(ptgg,wgg*0.84*wHalf->GetBinContent(theBin)); // gg: +/-16%
	pth2->Fill(ptgg,wgg*1.16*wHalf->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*0.84*wHalf->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*1.16*wHalf->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < 140. && mVBF > 105. && ptgg < 500.) {
	pth->Fill(ptVBF,wVBF);
        pth1->Fill(ptVBF,wVBF*1.02);             // VBF: +/-2%
	pth2->Fill(ptVBF,wVBF*0.98);
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
	pth1->Fill(ptgg,wgg*wHalf->GetBinContent(theBin));
        nloh1->Fill(nlogg,wgg*wHalf->GetBinContent(theBin));
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
        nloh2->Fill(nlozz,wzz);
      }
    } 

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < 140. && mcr > 105. && ptcr < 500.) {
	pth->Fill(ptcr);
        nloh->Fill(nlocr);
      }
    } 

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

    can.SaveAs("pt_systUnbRegion.gif");
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

    can.SaveAs("pt_systSingleZ.gif");
    return;
  }

  pth->Scale(1./pth->Integral());
  pth1->Scale(1./pth1->Integral());
  pth2->Scale(1./pth2->Integral());
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
  pth2->GetXaxis()->SetLabelColor(kWhite);
  pth2->GetYaxis()->SetTitle("Normalized entries");
  pth2->Draw();
  pth1->Draw("SAME");
  pth->Draw("SAME");
  
  TH1F* diffpt = (TH1F*)pth->Clone();
  TH1F* pullpt = (TH1F*)pth->Clone();
  diffpt->Add(pth2,pth,1,-1);
  pullpt->Divide(diffpt,pth);
  can.cd(2);
  gPad->SetTopMargin(0.0);
  pullpt->SetMinimum(-1.);
  pullpt->SetMaximum(1.);
  pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  pullpt->GetYaxis()->SetTitle("Relative difference");
  pullpt->Draw("E");
  
  sprintf(nameFile,"pt_syst%s.gif",nameSyst);
  can.SaveAs(nameFile);
  
  can.cd(1);
  gPad->SetBottomMargin(0.0);
  nloh1->GetXaxis()->SetLabelColor(kWhite);
  nloh1->GetYaxis()->SetTitle("Normalized entries");
  nloh2->GetXaxis()->SetLabelColor(kWhite);
  nloh2->GetYaxis()->SetTitle("Normalized entries");
  if (whichtype == 2) nloh2->Draw();
  else {
    nloh1->Draw();
    nloh2->Draw("SAME");
  }
  nloh->Draw("SAME");
  
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
  
  sprintf(nameFile,"nlo_syst%s.gif",nameSyst);
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
       
    sprintf(nameFile,"pt_balanced%s.gif",nameSyst);
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

    sprintf(nameFile,"nlo_balanced%s.gif",nameSyst);
    can.SaveAs(nameFile);
  }

  return;
}

