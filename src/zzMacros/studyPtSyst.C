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

static const int massRanges = 4;
static const int melaRanges = 3;

float ptVar(float pt, float m, bool overM = false) {
  if (overM) return pt/m;
  return pt;	
}

bool notVBFtagged(float njets) {
  if (njets >= 2) return false;
  return true;
}

void evalBinMigration(TH1F* def, TH1F* var, char fileName[200], bool syst = true, int binlimit = 20) { // i.e. 50 GeV
  
  char outpStr[200];
  sprintf(outpStr,"MCstatistics");
  // TH1F* def = (TH1F*)theFile->Get("ptH_Default");
  if (!def) return;
  float theErr;
  if (syst) {
    theErr = (def->Integral(binlimit,160)/def->Integral() - var->Integral(binlimit,160)/var->Integral())*(def->Integral()/def->Integral(binlimit,160));
    sprintf(outpStr,"%s",fileName);
  } else {
    theErr = 1./sqrt(def->GetEntries()*def->Integral(binlimit,160)/def->Integral());
  }
  cout << outpStr << " : " << int(theErr*10000)/100. << "%" << endl;
  return;
}

void adjustHistogram(TH1F* hist) {
  
  float a, b;
  for (Int_t iBin = 1; iBin <= hist->GetNbinsX() ; ++iBin) {
    if (hist->GetBinContent(iBin) < hist->GetBinError(iBin)) hist->SetBinError(iBin,0.9*hist->GetBinContent(iBin));
    if (iBin == 3) a = hist->GetBinContent(iBin);
    if (iBin == 4) b = hist->GetBinContent(iBin);
  }
}

void studyPtSyst(int mass = 125, int whichtype = 1, bool overM = false, 
		 bool also7TeV = true, bool moreFiles = false) { 
                                       // NEVER use moreFiles = true!

  // -5 - PDF : VBF
  // -4 - scales : VBF
  // -3 - effect of finite top mass
  // -2 - mu_Q in gg signal shape
  // -1 - fraction of VBF (NOT for VBF fraction fitting)
  // 0 - DEFAULTS: ALWAYS RUN FIRST
  // 1 - Z+X in background
  // 2 - Cross-check Z+X vs. Z+jets
  // 3 - low pT shape in background: nonblinded
  // 4 - low pT shape in background: single Z
  // 5 - adding/not adding ggZZ
  // 6 - PDF : ZZ
  // 7 - scales : ZZ

  bool withNLOMela = false;
  char nameFile[200] = "SelectedTree";
  char nameFile2[200] = "";
  if (withNLOMela) sprintf(nameFile,"angles");

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
 
  // READ TREES
  TChain* ggTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H%d.root",mass);  ggTree->Add(nameFile2);
    if (moreFiles) {
      /* ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H115.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H115.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H115.root"); 
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H120.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H120.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H120.root"); 
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H130.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H130.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H130.root"); 
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_H135.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_H135.root");
      ggTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_H135.root"); */
    } 
  }

  TChain* VBFTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VBFH%d.root",mass);    VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VBFH%d.root",mass);    VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VBFH%d.root",mass);    VBFTree->Add(nameFile2);
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH%d.root",mass);  VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH%d.root",mass);  VBFTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH%d.root",mass);  VBFTree->Add(nameFile2);
    if (moreFiles) {
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH126.root");
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH126.root");
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH126.root"); 
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VBFH127.root");
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VBFH127.root");
      VBFTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VBFH127.root"); 
    }
  }

  TChain* VHTree = new TChain(nameFile);
  if (also7TeV) {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VH%d.root",mass);   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VH%d.root",mass);   VHTree->Add(nameFile2);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VH%d.root",mass);   VHTree->Add(nameFile2);
    /* sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_VH130.root");
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_VH130.root");
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_VH130.root"); */
  } else {
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_VH%d.root",mass);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_VH%d.root",mass);
    sprintf(nameFile2,"root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_VH%d.root",mass); 
  }

  TChain* zzTree = new TChain(nameFile);
  if (also7TeV) {
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo4e.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo4mu.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2e2mu.root"); 
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ZZTo4tau.root");
  } else {
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4e.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4mu.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2mu.root"); 
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
    zzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ZZTo4tau.root");
  }

  TChain* ggzzTree = new TChain(nameFile);
  if (also7TeV) {
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4e/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/4mu/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  } else {
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4e/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/4mu/HZZ4lTree_ggZZ4l.root");
    ggzzTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  }
   
  TChain* dataTree = new TChain(nameFile);
  if (also7TeV) {
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleEle.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleMu.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/data/HZZ4lTree_DoubleOr.root");
  } else {
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleEle.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleMu.root");
    dataTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/data/HZZ4lTree_DoubleOr.root"); 
  }

  TChain* crTree = new TChain(nameFile);
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEMMssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMEEssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMMMssTree.root");
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");	
  crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMMMssTree.root"); 
  if (also7TeV) {
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEEEssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEMMssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMEEssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMMMssTree.root");
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEMMssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMEEssTree.root");	
    crTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMMMssTree.root");
  }

  TChain* crosTree = new TChain(nameFile);
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CREEMMosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMEEosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleEle_CRMMMMosTree.root");
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CREEMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");	
  crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DoubleOr_CRMMMMosTree.root");	
  if (also7TeV) {
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEEEosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CREEMMosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMEEosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleEle_CRMMMMosTree.root");
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CREEMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleMu_CRMMMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CREEMMosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMEEosTree.root");	
    crosTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR/CR/HZZ4lTree_DoubleOr_CRMMMMosTree.root");	
  }

  TChain* crzjTree = new TChain(nameFile);
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEEEssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CREEMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMEEssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50NoB_CRMMMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10NoB_CRMMMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M10B_CRMMMMssTree.root");
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CREEMMssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMEEssTree.root");	
  crzjTree->Add("root://lxcms02//data/Higgs/rootuplesOut/241112/PRODFSR_8TeV/CR/HZZ4lTree_DYJetsToLLTuneZ2M50B_CRMMMMssTree.root");

  float mgg, mVBF, mzz, mggzz, mdata, mcr, mcros, mzj, mvh;
  float wgg, wVBF, wzz, wggzz, wzj, wvh;
  float ptgg, ptVBF, ptzz, ptggzz, ptdata, ptcr, ptzj, ptvh;
  int njgg, njVBF, njzz, njggzz, njdata, njcr, njzj, njvh;
  float nlogg = -999.;
  float nloVBF = -999.;
  float nlozz = -999.; 
  float nloggzz = -999.;
  float nlodata = -999.; 
  float nlocr = -999.; 
  float nlozj = -999.; 
  float nlovh = -999.; 
  float genptgg;
  int gprIdvh; 

  ggTree->SetBranchAddress("ZZMass",&mgg);
  ggTree->SetBranchAddress("MC_weight",&wgg);
  ggTree->SetBranchAddress("ZZPt",&ptgg);
  ggTree->SetBranchAddress("NJets",&njgg);
  ggTree->SetBranchAddress("GenHPt",&genptgg);
  if (withNLOMela) ggTree->SetBranchAddress("melaLDWithPtY",&nlogg);
  else ggTree->SetBranchAddress("ZZLD",&nlogg);

  VBFTree->SetBranchAddress("ZZMass",&mVBF);
  VBFTree->SetBranchAddress("MC_weight",&wVBF);
  VBFTree->SetBranchAddress("ZZPt",&ptVBF);
  VBFTree->SetBranchAddress("NJets",&njVBF);
  if (withNLOMela) VBFTree->SetBranchAddress("melaLDWithPtY",&nloVBF);
  else VBFTree->SetBranchAddress("ZZLD",&nloVBF);  

  VHTree->SetBranchAddress("ZZMass",&mvh);
  VHTree->SetBranchAddress("MC_weight",&wvh);
  VHTree->SetBranchAddress("ZZPt",&ptvh);
  VHTree->SetBranchAddress("NJets",&njvh);
  VHTree->SetBranchAddress("genProcessId",&gprIdvh);
  if (withNLOMela) VHTree->SetBranchAddress("melaLDWithPtY",&nlovh);
  else VHTree->SetBranchAddress("ZZLD",&nlovh);

  zzTree->SetBranchAddress("ZZMass",&mzz);
  zzTree->SetBranchAddress("MC_weight",&wzz);
  zzTree->SetBranchAddress("ZZPt",&ptzz);
  zzTree->SetBranchAddress("NJets",&njzz);
  if (withNLOMela) zzTree->SetBranchAddress("melaLDWithPtY",&nlozz);
  else zzTree->SetBranchAddress("ZZLD",&nlozz);

  ggzzTree->SetBranchAddress("ZZMass",&mggzz);
  ggzzTree->SetBranchAddress("MC_weight",&wggzz);
  ggzzTree->SetBranchAddress("ZZPt",&ptggzz);
  ggzzTree->SetBranchAddress("NJets",&njggzz);
  if (withNLOMela) ggzzTree->SetBranchAddress("melaLDWithPtY",&nloggzz);
  else ggzzTree->SetBranchAddress("ZZLD",&nloggzz);

  dataTree->SetBranchAddress("ZZMass",&mdata);
  dataTree->SetBranchAddress("ZZPt",&ptdata);
  dataTree->SetBranchAddress("NJets",&njdata);
  if (withNLOMela) dataTree->SetBranchAddress("melaLDWithPtY",&nlodata);
  else dataTree->SetBranchAddress("ZZLD",&nlodata);

  crTree->SetBranchAddress("ZZMass",&mcr);
  crTree->SetBranchAddress("ZZPt",&ptcr);
  crTree->SetBranchAddress("NJets",&njcr);
  if (withNLOMela) crTree->SetBranchAddress("melaLDWithPtY",&nlocr);
  else crTree->SetBranchAddress("ZZLD",&nlocr);

  crosTree->SetBranchAddress("ZZMass",&mcros);

  crzjTree->SetBranchAddress("ZZMass",&mzj);
  crzjTree->SetBranchAddress("MC_weight",&wzj);
  crzjTree->SetBranchAddress("ZZPt",&ptzj);
  crzjTree->SetBranchAddress("NJets",&njzj);
  if (withNLOMela) crzjTree->SetBranchAddress("melaLDWithPtY",&nlozj);
  else crzjTree->SetBranchAddress("ZZLD",&nlozj);

  char nameSyst[200]; 
  char UcasePt[8] = "PT";
  if (overM) sprintf(UcasePt,"PToverM");
  char LcasePt[8] = "pt";
  if (overM) sprintf(LcasePt,"ptoverm");

  // cout << LcasePt << endl;

  sprintf(nameFile,"../../weightHisto_125GeV_8TeV_Default.root");
  if (also7TeV) sprintf(nameFile,"../../weightHisto_125GeV_7TeV_Default.root");
  TFile weightsSig(nameFile);
  TH1F* wHalf = (TH1F*)weightsSig.Get("wH");
  TH1F* wRenorm;
  if (also7TeV) {
    TFile weights8("../../weightHisto_125GeV_8TeV_Default.root");
    wRenorm = (TH1F*)weights8.Get("wH");
  } else {
    wRenorm = (TH1F*)weightsSig.Get("wH");
  }
  TFile weightsSig1("../../weightHisto_125GeV_8TeV_up.root");
  TH1F* wOne = (TH1F*)weightsSig1.Get("wH");
  TFile weightsSig2("../../weightHisto_125GeV_8TeV_down.root");
  TH1F* wQuar = (TH1F*)weightsSig2.Get("wH");

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
  if (overM == true) theMax = 3.2;

  float massLimits[massRanges+1] = {0.,0.,0.,0.,0.};   // fill after  
  float melaLimits[melaRanges+1] = {0.0,0.3,0.6,1.0}; 

  // Decide the width (un po' a vacca)
  if (mass <= 250) {massLimits[0] = mass-20.;  massLimits[massRanges] = mass+15.;}
  else if (mass <= 500) {massLimits[0] = mass*0.85;  massLimits[massRanges] = mass*1.15;}
  else if (mass <= 850) {massLimits[0] = mass*0.657;  massLimits[massRanges] = mass*1.571;}
  else {massLimits[0] = mass*0.4;  massLimits[massRanges] = mass*1.6;}
  for (int i = 1; i < massRanges; i++) {
    massLimits[i] = massLimits[0]+float(i/massRanges)*(massLimits[massRanges]-massLimits[0]);
  }

  // gg 
  TH1F* pth = new TH1F("pth","pt",nbins2,0.,theMax);
  pth->Sumw2();

  TH1F* pthMass[massRanges];
  TH1F* pthMela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"pthMass%d",i);
    pthMass[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    pthMass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"pthMela%d",i);
    pthMela[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    pthMela[i]->Sumw2();
  }

  // vbf
  TH1F* ptvbf = new TH1F("ptvbf","pt",nbins2,0.,theMax);
  ptvbf->Sumw2();
  
  TH1F* ptvbfMass[massRanges];
  TH1F* ptvbfMela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"%ptvbfMass%d",i);
    ptvbfMass[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    ptvbfMass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"ptvbfMela%d",i);
    ptvbfMela[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    ptvbfMela[i]->Sumw2();
  }

  // zz
  TH1F* pth1 = new TH1F("pth1","pt",nbins2,0.,theMax);
  pth1->Sumw2();
  TH1F* pth1Mass[massRanges];
  TH1F* pth1Mela[melaRanges];
  for (int i = 0; i < massRanges; i++) {
    sprintf(nameFile,"pth1Mass%d",i);
    pth1Mass[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    pth1Mass[i]->Sumw2();
  }
  for (int i = 0; i < melaRanges; i++) {
    sprintf(nameFile,"pth1Mela%d",i);
    pth1Mela[i] = new TH1F(nameFile,"pt",nbins2,0.,theMax);
    pth1Mela[i]->Sumw2();
  }

  // other
  TH1F* pth2 = new TH1F("pth2","pt",nbins2,0.,theMax);
  TH1F* pth3 = new TH1F("pth3","pt",nbins2,0.,theMax);
  TH1F* pth4 = new TH1F("pth4","pt",nbins2,0.,theMax);
  TH1F* pth5 = new TH1F("pth5","pt",nbins2,0.,theMax);
  TH1F* pth6 = new TH1F("pth6","pt",nbins2,0.,theMax);
  pth2->Sumw2();
  pth3->Sumw2();
  pth4->Sumw2();
  pth5->Sumw2();
  pth6->Sumw2();
 
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
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	int theBin = pdfh[0]->FindBin(ptVBF);
 	pth->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	pth1->Fill(ptVar(ptVBF,mVBF,overM),wVBF*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(ptVar(ptVBF,mVBF,overM),wVBF*ratiopdf[0]->GetBinContent(theBin));
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*ratiopdf[1]->GetBinContent(theBin));
	nloh2->Fill(nloVBF,wVBF*ratiopdf[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVar(ptVBF,mVBF,overM),wVBF*thisDiff); 
      }
    } 

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
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	int theBin = scaleh[0]->FindBin(ptVBF);
 	pth->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	pth1->Fill(ptVar(ptVBF,mVBF,overM),wVBF*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(ptVar(ptVBF,mVBF,overM),wVBF*ratioscale[0]->GetBinContent(theBin));
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*ratioscale[1]->GetBinContent(theBin));
	nloh2->Fill(nloVBF,wVBF*ratioscale[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVar(ptVBF,mVBF,overM),wVBF*thisDiff); 
      }
    } 

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

    // if (!also7TeV) can.SaveAs("pt_systTopMass_corr.gif");
    // if (!also7TeV) can.SaveAs("pt_systTopMass_corr.pdf");
    sprintf(nameSyst,"TopMass");

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	int theBin = wHalf->FindBin(genptgg);
        float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptVar(ptgg,mgg,overM),wgg*wHalf->GetBinContent(theBin));
	pth2->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
	ptvbf->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*wHalf->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
      }
    } 

    /* for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	pth2->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        nloh->Fill(nloVBF,wVBF);
	nloh2->Fill(nloVBF,wVBF);
      }
      }*/

    // return;

  }  

  if (whichtype == -2) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	int theBin = wHalf->FindBin(genptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
	pth1->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wQuar->GetBinContent(theBin)*wHalf->GetBinContent(theBin)/wRenorm->GetBinContent(theBin));
	pth2->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wOne->GetBinContent(theBin)*wHalf->GetBinContent(theBin)/wRenorm->GetBinContent(theBin));
	ptvbf->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wOne->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*newW*wOne->GetBinContent(theBin)*wHalf->GetBinContent(theBin)/wRenorm->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*wQuar->GetBinContent(theBin)*wHalf->GetBinContent(theBin)/wRenorm->GetBinContent(theBin));
      }
    } 

    sprintf(nameSyst,"Resummation");

    /* for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        pth1->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	pth2->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	ptvbf->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF);
	nloh2->Fill(nloVBF,wVBF);
      }
      }*/

  }

  if (whichtype == -1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	int theBin = wHalf->FindBin(genptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
	pth1->Fill(ptVar(ptgg,mgg,overM),wgg*newW*0.84*wHalf->GetBinContent(theBin));  // gg: +/-16%
	pth2->Fill(ptVar(ptgg,mgg,overM),wgg*newW*1.16*wHalf->GetBinContent(theBin));
	ptvbf->Fill(ptVar(ptgg,mgg,overM),wgg*newW*1.16*wHalf->GetBinContent(theBin));  
        nloh->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
	nloh1->Fill(nlogg,wgg*newW*0.84*wHalf->GetBinContent(theBin));
	nloh2->Fill(nlogg,wgg*newW*1.16*wHalf->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        pth1->Fill(ptVar(ptVBF,mVBF,overM),wVBF*1.02);             // VBF: +/-2%
	pth2->Fill(ptVar(ptVBF,mVBF,overM),wVBF*0.98);
	ptvbf->Fill(ptVar(ptVBF,mVBF,overM),wVBF*0.98);
        nloh->Fill(nloVBF,wVBF);
	nloh1->Fill(nloVBF,wVBF*1.02);
	nloh2->Fill(nloVBF,wVBF*0.98);
      }
    }

    sprintf(nameSyst,"VBFfraction");
 
  }

  if (whichtype == 0) {

    TH2F* ptkdgg = new TH2F("ptkdgg","pt vs mela",16,0.,1.,50.,0.,theMax);
    ptkdgg->Sumw2();
    TH2F* ptmhgg = new TH2F("ptmhgg","pt vs m",16,massLimits[0],massLimits[massRanges],50.,0.,theMax);
    ptmhgg->Sumw2();
    TH2F* ptmhggLarge = new TH2F("ptmhggLarge","pt vs m",90,100.,1000.,50.,0.,theMax);
    ptmhggLarge->Sumw2();

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	int theBin = wHalf->FindBin(genptgg);
	float newW = 1.;
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
        nloh->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
        ptkdgg->Fill(nlogg,ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin)); 
	for (int i = 0; i < melaRanges; i++) {
	  if (nlogg > melaLimits[i] && nlogg < melaLimits[i+1]) 
	    pthMela[i]->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mgg > massLimits[i] && mgg < massLimits[i+1] && ptgg < 400. && notVBFtagged(njgg)) 
	  pthMass[i]->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
      }      
      ptmhgg->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
      ptmhggLarge->Fill(mgg,ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
    } 

    sprintf(nameFile,"../../selRootFiles/%s_gg%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_gg%d_SEL_7TeV.root",UcasePt,mass);
    TFile f(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth->SetName(nameFile);
    pth->Scale(1./pth->Integral());
    adjustHistogram(pth);   pth->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      pthMass[i]->SetName(nameFile);
      pthMass[i]->Scale(1./pthMass[i]->Integral());
      adjustHistogram(pthMass[i]);   pthMass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      pthMela[i]->SetName(nameFile);
      pthMela[i]->Scale(1./pthMela[i]->Integral());
      adjustHistogram(pthMela[i]);   pthMela[i]->Write();
    }
    
    ptmhgg->SetName("ptVsM");   ptmhgg->Write();
    TProfile* tpgg = (TProfile*)ptmhgg->ProfileX();
    tpgg->SetName("ptVsMProf");   tpgg->Write();

    ptkdgg->SetName("ptVsMELA");   ptkdgg->Write();
    TProfile* tpkdgg = (TProfile*)ptkdgg->ProfileX();
    tpkdgg->SetName("ptVsMELAProf");   tpkdgg->Write();

    ptmhggLarge->SetName("ptVsMLarge");   ptmhggLarge->Write();
    TProfile* tpggLarge = (TProfile*)ptmhggLarge->ProfileX();
    tpggLarge->SetName("ptVsMLargeProf");   tpggLarge->Write();

    evalBinMigration(pth,pth,"",false);
    f.Close();

    TH2F* ptkdvbf = new TH2F("ptkdvbf","pt vs mela",16,0.,1.,50.,0.,theMax);
    ptkdvbf->Sumw2();
    TH2F* ptmhvbf = new TH2F("ptmhvbf","pt vs m",16,massLimits[0],massLimits[massRanges],50.,0.,theMax);
    ptmhvbf->Sumw2();
    TH2F* ptmhvbfLarge = new TH2F("ptmhvbfLarge","pt vs m",90,100.,1000.,50.,0.,theMax);
    ptmhvbfLarge->Sumw2();

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	ptvbf->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        nloh->Fill(nloVBF,wVBF);
        ptkdvbf->Fill(nloVBF,ptVar(ptVBF,mVBF,overM),wVBF);
	for (int i = 0; i < melaRanges; i++) {
	if (nloVBF > melaLimits[i] && nloVBF < melaLimits[i+1]) 
	  ptvbfMela[i]->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mVBF > massLimits[i] && mVBF < massLimits[i+1] && ptVBF < 400. && notVBFtagged(njVBF)) 
	  ptvbfMass[i]->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
      }
      ptmhvbf->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
      ptmhvbfLarge->Fill(mVBF,ptVar(ptVBF,mVBF,overM),wVBF);
    }

    sprintf(nameFile,"../../selRootFiles/%s_vbf%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_vbf%d_SEL_7TeV.root",UcasePt,mass);
    TFile f1(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    ptvbf->SetName(nameFile);
    ptvbf->Scale(1./ptvbf->Integral());
    adjustHistogram(ptvbf);   ptvbf->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      ptvbfMass[i]->SetName(nameFile);
      ptvbfMass[i]->Scale(1./ptvbfMass[i]->Integral());
      adjustHistogram(ptvbfMass[i]);   ptvbfMass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      ptvbfMela[i]->SetName(nameFile);
      ptvbfMela[i]->Scale(1./ptvbfMela[i]->Integral());
      adjustHistogram(ptvbfMela[i]);   ptvbfMela[i]->Write();
    }
    
    ptmhvbf->SetName("ptVsM");   ptmhvbf->Write();
    TProfile* tpvbf = (TProfile*)ptmhvbf->ProfileX();
    tpvbf->SetName("ptVsMProf");   tpvbf->Write();

    ptkdvbf->SetName("ptVsMELA");   ptkdvbf->Write();
    TProfile* tpkdvbf = (TProfile*)ptkdvbf->ProfileX();
    tpkdvbf->SetName("ptVsMELAProf");   tpkdvbf->Write();

    ptmhvbfLarge->SetName("ptVsMLarge");   ptmhvbfLarge->Write();
    TProfile* tpvbfLarge = (TProfile*)ptmhvbfLarge->ProfileX();
    tpvbfLarge->SetName("ptVsMLargeProf");   tpvbfLarge->Write();

    evalBinMigration(ptvbf,ptvbf,"",false);
    f1.Close();

    TH2F* ptkdzz = new TH2F("ptkdzz","pt vs mela",16,0.,1.,50.,0.,theMax);
    ptkdzz->Sumw2();
    TH2F* ptmhzz = new TH2F("ptmhzz","pt vs m",16,massLimits[0],massLimits[massRanges],50.,0.,theMax);
    ptmhzz->Sumw2();
    TH2F* ptmhzzLarge = new TH2F("ptmhzzLarge","pt vs m",90,100.,1000.,50.,0.,theMax);
    ptmhzzLarge->Sumw2();

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth1->Fill(ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz);
        nloh->Fill(nlozz,wzz);
        ptkdzz->Fill(nlozz,ptVar(ptzz,mzz,overM),wzz);
	for (int i = 0; i < melaRanges; i++) {
	  if (nlozz > melaLimits[i] && nlozz < melaLimits[i+1]) 
	    pth1Mela[i]->Fill(ptVar(ptzz,mzz,overM),wzz);
	}
      }
      for (int i = 0; i < massRanges; i++) {
	if (mzz > massLimits[i]-7. && mzz < massLimits[i+1]+7. && ptzz < 400. && notVBFtagged(njzz)) 
	  pth1Mass[i]->Fill(ptVar(ptzz,mzz,overM),wzz);
      }
      ptmhzz->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
      ptmhzzLarge->Fill(mzz,ptVar(ptzz,mzz,overM),wzz);
    } 

    sprintf(nameFile,"../../selRootFiles/%s_zz%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_zz%d_SEL_7TeV.root",UcasePt,mass);
    TFile f2(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth1->SetName(nameFile);
    pth1->Scale(1./pth1->Integral());
    adjustHistogram(pth1);   pth1->Write();
    for (int i = 0; i < massRanges; i++) {
      sprintf(nameFile,"%sH_Mass%d-%d",LcasePt,int(massLimits[i]),int(massLimits[i+1]));
      pth1Mass[i]->SetName(nameFile);
      pth1Mass[i]->Scale(1./pth1Mass[i]->Integral());
      adjustHistogram(pth1Mass[i]);   pth1Mass[i]->Write();
    }
    for (int i = 0; i < melaRanges; i++) {
      if (i == melaRanges - 1) sprintf(nameFile,"%sH_Mela0%d-10",LcasePt,int(melaLimits[i]*10));
      else sprintf(nameFile,"%sH_Mela0%d-0%d",LcasePt,int(melaLimits[i]*10),int(melaLimits[i+1]*10));
      pth1Mela[i]->SetName(nameFile);
      pth1Mela[i]->Scale(1./pth1Mela[i]->Integral());
      adjustHistogram(pth1Mela[i]);   pth1Mela[i]->Write();
    }
    
    ptmhzz->SetName("ptVsM");   ptmhzz->Write();
    TProfile* tpzz = (TProfile*)ptmhzz->ProfileX();
    tpzz->SetName("ptVsMProf");   tpzz->Write();

    ptkdzz->SetName("ptVsMELA");   ptkdzz->Write();
    TProfile* tpkdzz = (TProfile*)ptkdzz->ProfileX();
    tpkdzz->SetName("ptVsMELAProf");   tpkdzz->Write();

    ptmhzzLarge->SetName("ptVsMLarge");   ptmhzzLarge->Write();
    TProfile* tpzzLarge = (TProfile*)ptmhzzLarge->ProfileX();
    tpzzLarge->SetName("ptVsMLargeProf");   tpzzLarge->Write();

    evalBinMigration(pth1,pth1,"",false);
    f2.Close();

    TH2F* ptkdcr = new TH2F("ptkdcr","pt vs mela",16,0.,1.,50.,0.,theMax);
    ptkdcr->Sumw2();
    TH2F* ptmhcr = new TH2F("ptmhcr","pt vs m",16,massLimits[0],massLimits[massRanges],50.,0.,theMax);
    ptmhcr->Sumw2();
    TH2F* ptmhcrLarge = new TH2F("ptmhcrLarge","pt vs m",90,100.,1000.,50.,0.,theMax);
    ptmhcrLarge->Sumw2();

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < massLimits[massRanges] && mcr > massLimits[0] && ptcr < 400. && notVBFtagged(njcr)) {
	pth2->Fill(ptVar(ptcr,mcr,overM));
        nloh->Fill(nlocr);
        ptkdcr->Fill(nlocr,ptVar(ptcr,mcr,overM));
      }
      ptmhcr->Fill(mcr,ptVar(ptcr,mcr,overM));
      ptmhcrLarge->Fill(mcr,ptVar(ptcr,mcr,overM));
    } 

    sprintf(nameFile,"../../selRootFiles/%s_zx%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_zx%d_SEL_7TeV.root",UcasePt,mass);
    TFile f3(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth2->SetName(nameFile);
    pth2->Scale(1./pth2->Integral());
    adjustHistogram(pth2);   pth2->Write();

    ptmhcr->SetName("ptVsM");   ptmhcr->Write();
    TProfile* tpcr = (TProfile*)ptmhcr->ProfileX();
    tpcr->SetName("ptVsMProf");   tpcr->Write();

    ptkdcr->SetName("ptVsMELA");   ptkdcr->Write();
    TProfile* tpkdcr = (TProfile*)ptkdcr->ProfileX();
    tpkdcr->SetName("ptVsMELAProf");   tpkdcr->Write();

    ptmhcrLarge->SetName("ptVsMLarge");   ptmhcrLarge->Write();
    TProfile* tpcrLarge = (TProfile*)ptmhcrLarge->ProfileX();
    tpcrLarge->SetName("ptVsMLargeProf");   tpcrLarge->Write();

    evalBinMigration(pth2,pth2,"",false);
    f3.Close();

    for (Int_t iEvt = 0; iEvt < ggzzTree->GetEntries() ; ++iEvt) {
      ggzzTree->GetEntry(iEvt);
      if (mggzz < massLimits[massRanges] && mggzz > massLimits[0] && ptggzz < 400. && notVBFtagged(njggzz)) {
	pth3->Fill(ptVar(ptggzz,mggzz,overM),wggzz);
        nloh->Fill(nloggzz,wggzz);
      }
    } 

    sprintf(nameFile,"../../selRootFiles/%s_ggzz%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_ggzz%d_SEL_7TeV.root",UcasePt,mass);
    TFile f4(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth3->SetName(nameFile);
    pth3->Scale(1./pth3->Integral());
    adjustHistogram(pth3);   pth3->Write();
    evalBinMigration(pth3,pth3,"",false);
    f4.Close();

    for (Int_t iEvt = 0; iEvt < VHTree->GetEntries() ; ++iEvt) {
      VHTree->GetEntry(iEvt);
      if (mvh < massLimits[massRanges] && mvh > massLimits[0] && ptvh < 400. && notVBFtagged(njvh)) {
	if (gprIdvh == 24) {
	  pth4->Fill(ptVar(ptvh,mvh,overM),wvh);
	  nloh->Fill(nlovh,wvh);
	} else if (gprIdvh == 26) {
	  pth5->Fill(ptVar(ptvh,mvh,overM),wvh);
	  nloh->Fill(nlovh,wvh);
	} else {
	  pth6->Fill(ptVar(ptvh,mvh,overM),wvh);
	  nloh->Fill(nlovh,wvh);
	}
      }
    } 

    sprintf(nameFile,"../../selRootFiles/%s_wh%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_wh%d_SEL_7TeV.root",UcasePt,mass);
    TFile f5(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth4->SetName(nameFile);
    pth4->Scale(1./pth4->Integral());
    adjustHistogram(pth4);   pth4->Write();
    evalBinMigration(pth4,pth4,"",false);
    f5.Close();

    sprintf(nameFile,"../../selRootFiles/%s_zh%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_zh%d_SEL_7TeV.root",UcasePt,mass);
    TFile f6(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth5->SetName(nameFile);
    pth5->Scale(1./pth5->Integral());
    adjustHistogram(pth5);   pth5->Write();
    evalBinMigration(pth5,pth5,"",false);
    f6.Close();
 
    sprintf(nameFile,"../../selRootFiles/%s_tth%d_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_tth%d_SEL_7TeV.root",UcasePt,mass);
    TFile f7(nameFile,"RECREATE");
    sprintf(nameFile,"%sH_Default",LcasePt);   
    pth6->SetName(nameFile);
    pth6->Scale(1./pth6->Integral());
    adjustHistogram(pth6);   pth6->Write();
    evalBinMigration(pth6,pth6,"",false);
    f7.Close();

  }

  if (whichtype == 1) {

    for (Int_t iEvt = 0; iEvt < ggTree->GetEntries() ; ++iEvt) {
      ggTree->GetEntry(iEvt);
      if (mgg < massLimits[massRanges] && mgg > massLimits[0] && ptgg < 400. && notVBFtagged(njgg)) {
	int theBin = wHalf->FindBin(genptgg);
	float newW = 1.;
	// if (ptgg > 90.) 
	newW = ratio->GetBinContent(ratio->FindBin(ptgg));
	pth1->Fill(ptVar(ptgg,mgg,overM),wgg*newW*wHalf->GetBinContent(theBin));
        nloh1->Fill(nlogg,wgg*newW*wHalf->GetBinContent(theBin));
      }
    } 

    for (Int_t iEvt = 0; iEvt < VBFTree->GetEntries() ; ++iEvt) {
      VBFTree->GetEntry(iEvt);
      if (mVBF < massLimits[massRanges] && mVBF > massLimits[0] && ptVBF < 400. && notVBFtagged(njVBF)) {
	pth1->Fill(ptVar(ptVBF,mVBF,overM),wVBF);
        nloh1->Fill(nloVBF,wVBF);
      }
    }

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth2->Fill(ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz);
        nloh2->Fill(nlozz,wzz);
      }
    } 

    /* sprintf(nameFile,"%s_zz_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"%s_zz_SEL_7TeV.root",UcasePt,mass);
    TFile f2(nameFile,"UPDATE");
    pth2->SetName("ptH_Default");
    pth2->Write();
    evalBinMigration(pth2,pth2,"",false);
    f2.Close(); */

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < massLimits[massRanges] && mcr > massLimits[0] && ptcr < 400. && notVBFtagged(njcr)) {
	pth->Fill(ptVar(ptcr,mcr,overM));
        nloh->Fill(nlocr);
      }
    } 

    /* sprintf(nameFile,"%s_zx_SEL_8TeV.root",UcasePt,mass);
    if (also7TeV) sprintf(nameFile,"%s_zx_SEL_7TeV.root",UcasePt,mass);
    TFile f3(nameFile,"UPDATE");
    pth->SetName("ptH_Default");
    pth->Write();
    evalBinMigration(pth,pth,"",false);
    f3.Close(); */

    sprintf(nameSyst,"ZplusX");
 
  }

  if (whichtype == 2) {

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (ptcr < 400. && notVBFtagged(njcr)) {          // do it in the whole mass range
	pth->Fill(ptVar(ptcr,mcr,overM));
        nloh->Fill(nlocr);
      }
    } 

    for (Int_t iEvt = 0; iEvt < crzjTree->GetEntries() ; ++iEvt) {
      crzjTree->GetEntry(iEvt);
      if (ptzj < 400. && notVBFtagged(njzj)) {
	pth2->Fill(ptVar(ptzj,mzj,overM),wzj);
	ptvbf->Fill(ptVar(ptzj,mzj,overM),wzj);
        nloh2->Fill(nlozj,wzj);	
      }
    }

    sprintf(nameSyst,"ZjetsXcheck");
 
  }

  if (whichtype == 3) {

    // DO NOT USE UGLY HISTOGRAMS, USE THE TREES INSTEAD
    /* sprintf(nameSyst,"histograms8TeV.root");
    if (also7TeV) sprintf(nameSyst,"histograms7TeV.root");
    TFile fileMike(nameSyst);
    fileMike.ls();

    sprintf(nameSyst,"DATA");
    if (!also7TeV) sprintf(nameSyst,"tmpTH1");
    TH1F* dataUnb = (TH1F*)fileMike.Get(nameSyst);
    dataUnb->GetXaxis()->SetRange(1,20);
    TH1F* zzUnb = (TH1F*)fileMike.Get("ZZ");
    zzUnb->Scale(13.6/zzUnb->GetBinContent(4));
    zzUnb->GetXaxis()->SetRange(1,20);
    TH1F* zxUnb = (TH1F*)fileMike.Get("FAKES");
    zxUnb->Scale(1.7/zxUnb->GetBinContent(4));
    zxUnb->GetXaxis()->SetRange(1,20); */

    float unbMaxBin = 200.;
    if (overM) unbMaxBin = 0.8;

    TH1F* dataUnb = new TH1F("dataUnb","data",20,0.,unbMaxBin);
    TH1F* zzUnb = new TH1F("zzUnb","zz",20,0.,unbMaxBin);
    TH1F* zxUnb = new TH1F("zxUnb","zx",20,0.,unbMaxBin);
    dataUnb->Sumw2();
    zzUnb->Sumw2();
    zxUnb->Sumw2();

    for (Int_t iEvt = 0; iEvt < dataTree->GetEntries() ; ++iEvt) {
      dataTree->GetEntry(iEvt);
      if (mdata < 400. && mdata > massLimits[massRanges] && ptdata < 200. && notVBFtagged(njdata)) {
	dataUnb->Fill(ptVar(ptdata,mdata,overM));
      }
    } 

    float lumi = 12.2;
    if (also7TeV) lumi = 5.1;

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < 400. && mzz > massLimits[massRanges] && ptzz < 200. && notVBFtagged(njzz)) {
	zzUnb->Fill(ptVar(ptzz,mzz,overM),wzz*lumi);
      }
    } 

    for (Int_t iEvt = 0; iEvt < crTree->GetEntries() ; ++iEvt) {
      crTree->GetEntry(iEvt);
      if (mcr < 400. && mcr > massLimits[massRanges] && ptcr < 200. && notVBFtagged(njcr)) {
	zxUnb->Fill(ptVar(ptcr,mcr,overM));
      }
    } 

    // Subtract fakes 
    float weightZX = (4.193/20.516)*(zzUnb->Integral()/zxUnb->Integral());
    dataUnb->Add(dataUnb,zxUnb,1,-weightZX);
    dataUnb->SetMarkerColor(1);    dataUnb->SetLineColor(1);    dataUnb->SetLineWidth(2);
    zzUnb->SetMarkerColor(2);   zzUnb->SetLineColor(2);   zzUnb->SetLineWidth(2);
    // zxUnb->SetMarkerColor(3);   zxUnb->SetLineColor(3);   zxUnb->SetLineWidth(2);

    can.cd(1);
    gPad->SetBottomMargin(0.0);
    dataUnb->GetXaxis()->SetLabelColor(kWhite);
    dataUnb->GetYaxis()->SetTitle("Entries");
    dataUnb->Draw();
    zzUnb->Draw("SAME");
    // zxUnb->Draw("SAME");

    TH1F* diffpt = (TH1F*)zxUnb->Clone();
    TH1F* pullpt = (TH1F*)zxUnb->Clone();
    diffpt->Add(dataUnb,zzUnb,1,-1);
    pullpt->Divide(diffpt,dataUnb);
    // remove outliers
    float xp[20];     float yp[20];
    float xperr[20];  float yperr[20];
    for (Int_t iBin = 1; iBin <= pullpt->GetNbinsX() ; ++iBin) {
      xperr[iBin-1] = 0.;
      if (fabs(pullpt->GetBinContent(iBin)) > 2.) {
	xp[iBin-1] = 999.;   yp[iBin-1] = 0.;     yperr[iBin-1] = 0.;
      } else {
	xp[iBin-1] = pullpt->GetBinCenter(iBin);   
	yp[iBin-1] = pullpt->GetBinContent(iBin);     
	yperr[iBin-1] = pullpt->GetBinError(iBin);;
      }
    }
    TGraphErrors* pullg = new TGraphErrors(20,xp,yp,xperr,yperr);

    can.cd(2);
    gPad->SetTopMargin(0.0);
    pullpt->SetMinimum(-2.);
    pullpt->SetMaximum(2.);
    pullpt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    pullpt->GetYaxis()->SetTitle("Relative difference");
    pullpt->Draw("E");
    TF1 *mypol1 = new TF1("mypol1","pol1");
    pullg->Fit("mypol1","","",0.,9.*unbMaxBin/20.);
    pullg->Draw("PSAME");

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(ptVar(ptzz,mzz,overM),wzz*(1. + mypol1->Eval(ptVar(ptzz,mzz,overM))));
      }
    } 
    
    sprintf(nameSyst,"UnbRegion");

    sprintf(nameFile,"%s_systUnbRegion.gif",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    sprintf(nameFile,"%s_systUnbRegion.pdf",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    // return;

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
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(ptVar(ptzz,mzz,overM),wzz*(1. + mypol1->Eval(ptzz)));
      }
    } 

    sprintf(nameSyst,"SingleZ");    
    
    sprintf(nameFile,"%s_systSingleZ.gif",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);
    sprintf(nameFile,"%s_systSingleZ.pdf",LcasePt);
    if (!also7TeV) can.SaveAs(nameFile);

  }

  if (whichtype == 5) {

    for (Int_t iEvt = 0; iEvt < zzTree->GetEntries() ; ++iEvt) {
      zzTree->GetEntry(iEvt);
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	pth->Fill(ptVar(ptzz,mzz,overM),wzz);
        nloh->Fill(nlozz,wzz);
	pth2->Fill(ptVar(ptzz,mzz,overM),wzz);
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz);
        nloh2->Fill(nlozz,wzz);
      }
    } 

    for (Int_t iEvt = 0; iEvt < ggzzTree->GetEntries() ; ++iEvt) {
      ggzzTree->GetEntry(iEvt);
      if (mggzz < massLimits[massRanges] && mggzz > massLimits[0] && ptggzz < 400. && notVBFtagged(njggzz)) {
	pth1->Fill(ptVar(ptggzz,mggzz,overM),wggzz);
        nloh1->Fill(nloggzz,wggzz);
	pth->Fill(ptVar(ptggzz,mggzz,overM),wggzz);
        nloh->Fill(nloggzz,wggzz);
      }
    } 

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
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	int theBin = pdfh[0]->FindBin(ptzz);
 	pth->Fill(ptVar(ptzz,mzz,overM),wzz);
	pth1->Fill(ptVar(ptzz,mzz,overM),wzz*ratiopdf[1]->GetBinContent(theBin));
	pth2->Fill(ptVar(ptzz,mzz,overM),wzz*ratiopdf[0]->GetBinContent(theBin));
        nloh->Fill(nlozz,wzz);
	nloh1->Fill(nlozz,wzz*ratiopdf[1]->GetBinContent(theBin));
	nloh2->Fill(nlozz,wzz*ratiopdf[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 2; i++) {
	  if (fabs(ratiopdf[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratiopdf[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz*thisDiff); 
      }
    } 

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
      if (mzz < massLimits[massRanges] && mzz > massLimits[0] && ptzz < 400. && notVBFtagged(njzz)) {
	int theBin = scaleh[0]->FindBin(ptzz);
 	pth->Fill(ptVar(ptzz,mzz,overM),wzz);
	pth1->Fill(ptVar(ptzz,mzz,overM),wzz*ratioscale[1]->GetBinContent(theBin));
	pth2->Fill(ptVar(ptzz,mzz,overM),wzz*ratioscale[0]->GetBinContent(theBin));
        nloh->Fill(nlozz,wzz);
	nloh1->Fill(nlozz,wzz*ratioscale[1]->GetBinContent(theBin));
	nloh2->Fill(nlozz,wzz*ratioscale[0]->GetBinContent(theBin));
	float thisDiff = 1.;
	for (int i = 0; i < 6; i++) {
	  if (fabs(ratioscale[i]->GetBinContent(theBin) - 1.) > fabs(thisDiff - 1.)) thisDiff = ratioscale[i]->GetBinContent(theBin);
	}
	ptvbf->Fill(ptVar(ptzz,mzz,overM),wzz*thisDiff); 
      }
    } 

    sprintf(nameSyst,"scale-ZZ");

  }

  if (whichtype == 0) return;

  pth->Scale(1./pth->Integral());
  pth1->Scale(1./pth1->Integral());
  pth2->Scale(1./pth2->Integral());
  ptvbf->Scale(1./ptvbf->Integral());
  if (withNLOMela) {
    nloh->Scale(1./nloh->Integral());
    nloh1->Scale(1./nloh1->Integral());
    nloh2->Scale(1./nloh2->Integral());
  }

  // Save histos
  
  char whichSample[4] = "gg";
  if (whichtype == -4 || whichtype == -5) sprintf(whichSample,"vbf");
  if (whichtype > 0) sprintf(whichSample,"zz");

  sprintf(nameFile,"../../selRootFiles/%s_%s%d_SEL_8TeV.root",UcasePt,whichSample,mass);
  if (also7TeV) sprintf(nameFile,"../../selRootFiles/%s_%s%d_SEL_7TeV.root",UcasePt,whichSample,mass);
  TFile f1(nameFile,"UPDATE");
  sprintf(nameFile,"%sH_%s",LcasePt,nameSyst);
  // cout << whichtype << " " << nameFile << endl;
  if (whichtype == -4 || whichtype == -5 || whichtype == 6 || whichtype == 7) {
    ptvbf->SetName(nameFile);
    adjustHistogram(ptvbf);   ptvbf->Write();
  } else if (whichtype == -3 || whichtype == 3 || whichtype == 4 || whichtype == 5) {
    pth->SetName(nameFile);
    adjustHistogram(pth);   pth->Write();
  } else if (whichtype == -2) {
    pth2->SetName(nameFile);
    adjustHistogram(pth2);  pth2->Write();
  }
  sprintf(nameFile,"%sH_Default",LcasePt);
  TH1F* def = (TH1F*)f1.Get(nameFile);
  if (whichtype == -4 || whichtype == -5 || whichtype == 6 || whichtype == 7) evalBinMigration(def,ptvbf,nameSyst,true);
  else if (whichtype == -3 || whichtype == 3 ||  whichtype == 4 || whichtype == 5) evalBinMigration(def,pth,nameSyst,true);
  else if (whichtype == -2) evalBinMigration(def,pth2,nameSyst,true);
  f1.Close();

  if (whichtype == 3 || whichtype == 4) return;

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
  
  sprintf(nameFile,"%s_syst%s.gif",LcasePt,nameSyst);
  if (!also7TeV || (whichtype == 5 && also7TeV) ) can.SaveAs(nameFile);
  sprintf(nameFile,"%s_syst%s.pdf",LcasePt,nameSyst);
  if (!also7TeV || (whichtype == 5 && also7TeV) ) can.SaveAs(nameFile);

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
  
  sprintf(nameFile,"nlo_syst%s.gif",nameSyst);
  if (!also7TeV || (whichtype == 5 && also7TeV)) can.SaveAs(nameFile);
  sprintf(nameFile,"nlo_syst%s.pdf",nameSyst);
  if (!also7TeV || (whichtype == 5 && also7TeV)) can.SaveAs(nameFile);

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
       
    sprintf(nameFile,"%s_balanced%s.pdf",LcasePt,nameSyst);
    if (!also7TeV) can.SaveAs(nameFile);
    
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
    if (!also7TeV) can.SaveAs(nameFile);
  }

  return;
}

