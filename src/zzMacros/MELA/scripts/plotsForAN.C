#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace RooFit;
using namespace std;

// = = = = = = = = = = = = = = = = = = = = 
// plots show 8 projections for sig & bkg
// = = = = = = = = = = = = = = = = = = = =

void anglesAndMasses(int index=0){

  double lowM[3] = {100,180,300};
  double highM[3]= {135,220,500};
  int mH[3]   = {120,200,400};

  char fileName[150];
  sprintf(fileName,"../datafiles/ZZ*AnalysisTree_H%i_withDiscriminants.root",mH[index]);
    

  TChain* sigTree =new TChain("angles");
  sigTree->Add(fileName);
  TChain* bkgTree =new TChain("angles");
  bkgTree->Add("../datafiles/ZZ*AnalysisTree_ZZTo*_withDiscriminants.root");

  if(!sigTree || sigTree->GetEntries()<=0 || !bkgTree || bkgTree->GetEntries()<=0) 
    return;

  TH1F* hsig_m1 = new TH1F("hsig_m1",";m_{Z1};",35,50,120);
  TH1F* hbkg_m1 = new TH1F("hbkg_m1",";m_{Z1};",35,50,120);
  hsig_m1->Sumw2();
  hbkg_m1->Sumw2();
  TH1F* hsig_m2 = new TH1F("hsig_m2",";m_{Z2};",55,10,120);
  TH1F* hbkg_m2 = new TH1F("hbkg_m2",";m_{Z2};",55,10,120);
  hsig_m2->Sumw2();
  hbkg_m2->Sumw2();
  TH1F* hsig_h2 = new TH1F("hsig_h2",";cos#theta_{2};",20,-1,1);
  TH1F* hbkg_h2 = new TH1F("hbkg_h2",";cos#theta_{2};",20,-1,1);
  hsig_h2->Sumw2();
  hbkg_h2->Sumw2();
  TH1F* hsig_h1 = new TH1F("hsig_h1",";cos#theta_{1};",20,-1,1);
  TH1F* hbkg_h1 = new TH1F("hbkg_h1",";cos#theta_{1};",20,-1,1);
  hsig_h1->Sumw2();
  hbkg_h1->Sumw2();
  TH1F* hsig_hs = new TH1F("hsig_hs",";cos#theta^{*};",20,-1,1);
  TH1F* hbkg_hs = new TH1F("hbkg_hs",";cos#theta^{*};",20,-1,1);
  hsig_hs->Sumw2();
  hbkg_hs->Sumw2();
  TH1F* hsig_phi = new TH1F("hsig_phi",";#Phi;",20,-3.1415,3.1415);
  TH1F* hbkg_phi = new TH1F("hbkg_phi",";#Phi;",20,-3.1415,3.1415);
  hsig_phi->Sumw2();
  hbkg_phi->Sumw2();
  TH1F* hsig_phi1 = new TH1F("hsig_phi1",";#Phi_{1};",20,-3.1415,3.1415);
  TH1F* hbkg_phi1 = new TH1F("hbkg_phi1",";#Phi_{1};",20,-3.1415,3.1415);
  hsig_phi1->Sumw2();
  hbkg_phi1->Sumw2();

  float mzz,h1,h2,hs,phi,phi1,m1,m2,w;

  sigTree->SetBranchAddress("zzmass",&mzz);
  sigTree->SetBranchAddress("z1mass",&m1);
  sigTree->SetBranchAddress("z2mass",&m2);
  sigTree->SetBranchAddress("costheta1",&h1);
  sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("phi",&phi);
  sigTree->SetBranchAddress("phistar1",&phi1);
  sigTree->SetBranchAddress("MC_weight",&w);

  bkgTree->SetBranchAddress("zzmass",&mzz);
  bkgTree->SetBranchAddress("z1mass",&m1);
  bkgTree->SetBranchAddress("z2mass",&m2);
  bkgTree->SetBranchAddress("costheta1",&h1);
  bkgTree->SetBranchAddress("costheta2",&h2);
  bkgTree->SetBranchAddress("costhetastar",&hs);
  bkgTree->SetBranchAddress("phi",&phi);
  bkgTree->SetBranchAddress("phistar1",&phi1);
  bkgTree->SetBranchAddress("MC_weight",&w);
  
  for(int i=0; i<sigTree->GetEntries(); i++){

    sigTree->GetEntry(i);
    
    //cout << mzz << " " << m2 << " " << m1 << " " << h1 << " " << h2 << " " << hs << endl;

    if(mzz>lowM[index]&&mzz<highM[index]){

      hsig_m1->Fill(m1,w);
      hsig_m2->Fill(m2,w);
      hsig_h1->Fill(h1,w);
      hsig_h2->Fill(h2,w);
      hsig_hs->Fill(hs,w);
      hsig_phi->Fill(phi,w);
      hsig_phi1->Fill(phi1,w);

    }

  }

  vector<TH1F*> histoSig;
  histoSig.push_back(hsig_m1);
  histoSig.push_back(hsig_m2);
  histoSig.push_back(hsig_h1);
  histoSig.push_back(hsig_h2);
  histoSig.push_back(hsig_hs);
  histoSig.push_back(hsig_phi);
  histoSig.push_back(hsig_phi1);
  
  for(int i=0; i<bkgTree->GetEntries(); i++){
    
    bkgTree->GetEntry(i);

      //cout << mzz << " " << m2 << " " << m1 << " " << h1 << " " << h2 << " " << hs << endl;
      
      if(mzz>lowM[index]&&mzz<highM[index]){
	
	hbkg_m1->Fill(m1,w);
	hbkg_m2->Fill(m2,w);
	hbkg_h1->Fill(h1,w);
	hbkg_h2->Fill(h2,w);
	hbkg_hs->Fill(hs,w);
	hbkg_phi->Fill(phi,w);
	hbkg_phi1->Fill(phi1,w);
	
      }

    }

  vector<TH1F*> histoBkg;
  histoBkg.push_back(hbkg_m1);
  histoBkg.push_back(hbkg_m2);
  histoBkg.push_back(hbkg_h1);
  histoBkg.push_back(hbkg_h2);
  histoBkg.push_back(hbkg_hs);
  histoBkg.push_back(hbkg_phi);
  histoBkg.push_back(hbkg_phi1);

  for(int i=0; i<histoBkg.size(); i++){

    histoSig[i]->Scale(1/histoSig[i]->Integral());
    histoSig[i]->GetYaxis()->SetRangeUser(0,1.3*(histoSig[i]->GetMaximum()+histoSig[i]->GetBinError(histoSig[i]->GetMaximumBin())));
    histoSig[i]->SetLineColor(2);
    histoSig[i]->SetLineWidth(2);
    
    histoBkg[i]->Scale(1/histoBkg[i]->Integral());
    histoBkg[i]->GetYaxis()->SetRangeUser(0,1.3*(histoSig[i]->GetMaximum()+histoSig[i]->GetBinError(histoSig[i]->GetMaximumBin())));
    histoBkg[i]->SetLineColor(4);
    histoBkg[i]->SetLineWidth(2);
    histoBkg[i]->SetLineStyle(2);

  }

  TCanvas *can_m1 = new TCanvas("can_m1","can_m1",400,400);
  TCanvas *can_m2 = new TCanvas("can_m2","can_m2",400,400);
  TCanvas *can_h1 = new TCanvas("can_h1","can_h1",400,400);
  TCanvas *can_h2 = new TCanvas("can_h2","can_h2",400,400);
  TCanvas *can_hs = new TCanvas("can_hs","can_hs",400,400);
  TCanvas *can_phi = new TCanvas("can_phi","can_phi",400,400);
  TCanvas *can_phi1 = new TCanvas("can_phi1","can_phi1",400,400);
  
  vector<TCanvas*> canvases;
  canvases.push_back(can_m1);
  canvases.push_back(can_m2);
  canvases.push_back(can_h1);
  canvases.push_back(can_h2);
  canvases.push_back(can_hs);
  canvases.push_back(can_phi);
  canvases.push_back(can_phi1);
  
  vector<string> names;
  names.push_back("m1");
  names.push_back("m2");
  names.push_back("h1");
  names.push_back("h2");
  names.push_back("hs");
  names.push_back("phi");
  names.push_back("phi1");

  char temp[50];

  for(int i=0; i<canvases.size(); i++){

    canvases[i]->cd();
    
    if(histoBkg[i]->GetMaximum()>histoSig[i]->GetMaximum()){
      histoBkg[i]->Draw("HISTE");
      histoSig[i]->Draw("SAMEHISTE");
    }else{
      histoSig[i]->Draw("HISTE");
      histoBkg[i]->Draw("SAMEHISTE");
    }
    
    sprintf(temp,"compareSignalvsBackground_H%i_%s.eps",mH[index],names[i].c_str());
    canvases[i]->SaveAs(temp);

    sprintf(temp,"compareSignalvsBackground_H%i_%s.png",mH[index],names[i].c_str());
    canvases[i]->SaveAs(temp);
    
  }

}

// = = = = = = = = = = = = = = = = = 
// plots 2D heat maps of templates
// = = = = = = = = = = = = = = = = = 

void MELAtemplate(char* channel="4mu",bool lowMass=true,bool plotSmooth=true){

  char fileName[150];
  sprintf(fileName,"../datafiles/Dsignal_%s.root",channel);
  TFile* sigFile = new TFile(fileName);
  sprintf(fileName,"../datafiles/Dsignal_PS_%s.root",channel);
  TFile* sigPsFile = new TFile(fileName);
  sprintf(fileName,"../datafiles/Dbackground_qqZZ_%s.root",channel);
  TFile* qqZZFile = new TFile(fileName);
  sprintf(fileName,"../datafiles/Dbackground_ggZZ_%s.root",channel);
  TFile* ggZZFile = new TFile(fileName);

  TH2F* sigTemplate ;
  TH2F* sigPsTemplate ;
  TH2F* qqZZTemplate;
  TH2F* ggZZTemplate;

  if(plotSmooth){
    sigTemplate  = (TH2F*) sigFile->Get("h_mzzD");
    sigPsTemplate  = (TH2F*) sigPsFile->Get("h_mzzD");
    qqZZTemplate = (TH2F*) qqZZFile->Get("h_mzzD");
    ggZZTemplate = (TH2F*) ggZZFile->Get("h_mzzD");
  }else{
    sigTemplate  = (TH2F*) sigFile->Get("oldTemp");
    sigPsTemplate  = (TH2F*) sigPsFile->Get("oldTemp");
    qqZZTemplate = (TH2F*) qqZZFile->Get("oldTemp");
    ggZZTemplate = (TH2F*) ggZZFile->Get("oldTemp");
  }
  if(lowMass){
    sigTemplate->GetXaxis()->SetRangeUser(100,180);
    sigPsTemplate->GetXaxis()->SetRangeUser(100,180);
    qqZZTemplate->GetXaxis()->SetRangeUser(100,180);
    ggZZTemplate->GetXaxis()->SetRangeUser(100,180);
  }else{
    sigTemplate->GetXaxis()->SetRangeUser(180,800);
    sigPsTemplate->GetXaxis()->SetRangeUser(180,800);
    qqZZTemplate->GetXaxis()->SetRangeUser(180,800); 
    ggZZTemplate->GetXaxis()->SetRangeUser(180,800);
  }
 
  sigTemplate->GetXaxis()->SetTitle("m_{4l}");
  sigTemplate->GetYaxis()->SetTitle("D");
  sigPsTemplate->GetXaxis()->SetTitle("m_{4l}");
  sigPsTemplate->GetYaxis()->SetTitle("D");
  qqZZTemplate->GetXaxis()->SetTitle("m_{4l}");
  qqZZTemplate->GetYaxis()->SetTitle("D");
  ggZZTemplate->GetXaxis()->SetTitle("m_{4l}");
  ggZZTemplate->GetYaxis()->SetTitle("D");
  

  TCanvas* canSig = new TCanvas("canSig","canSig",400,400);
  sigTemplate->Draw("COL");
  TCanvas* canSigPs = new TCanvas("canSigPs","canSigPs",400,400);
  sigPsTemplate->Draw("COL");
  TCanvas* canqqZZ = new TCanvas("canqqZZ","canqqZZ",400,400);
  qqZZTemplate->Draw("COL");
  TCanvas* canggZZ = new TCanvas("canggZZ","canggZZ",400,400);
  ggZZTemplate->Draw("COL");

  if(lowMass)
    sprintf(fileName,"MELAtemplate%s_signal_%s_lowMass.eps",(plotSmooth)?"Smooth":"",channel);
  else
    sprintf(fileName,"MELAtemplate%s_signal_%s_highMass.eps",(plotSmooth)?"Smooth":"",channel);

  canSig->SaveAs(fileName);

  if(lowMass)
    sprintf(fileName,"MELAtemplate%s_signal_PS_%s_lowMass.eps",(plotSmooth)?"Smooth":"",channel);
  else
    sprintf(fileName,"MELAtemplate%s_signal_PS_%s_highMass.eps",(plotSmooth)?"Smooth":"",channel);

  canSigPs->SaveAs(fileName);

  if(lowMass)
    sprintf(fileName,"MELAtemplate%s_qqZZbackground_%s_lowMass.eps",(plotSmooth)?"Smooth":"",channel);
  else
    sprintf(fileName,"MELAtemplate%s_qqZZbackground_%s_highMass.eps",(plotSmooth)?"Smooth":"",channel);

  canqqZZ->SaveAs(fileName);

  if(lowMass)
    sprintf(fileName,"MELAtemplate%s_ggZZbackground_%s_lowMass.eps",(plotSmooth)?"Smooth":"",channel);
  else
    sprintf(fileName,"MELAtemplate%s_ggZZbackground_%s_highMass.eps",(plotSmooth)?"Smooth":"",channel);

  canggZZ->SaveAs(fileName);

}

void makeAllMELAtemplate(){

  MELAtemplate("4mu",false,true);
  MELAtemplate("4mu",true,true);
  MELAtemplate("4mu",false,false);
  MELAtemplate("4mu",true,false);
  
  MELAtemplate("4e",false,true);
  MELAtemplate("4e",true,true);
  MELAtemplate("4e",false,false);
  MELAtemplate("4e",true,false);

  MELAtemplate("2e2mu",false,true);
  MELAtemplate("2e2mu",true,true);
  MELAtemplate("2e2mu",false,false);
  MELAtemplate("2e2mu",true,false);

}

// = = = = = = = = = = = = = = = = = = = = = =
// compare tempates to individual MC samples
// = = = = = = = = = = = = = = = = = = = = = =

void cocktailSyst2(int mass, char* channel, char* tempFileName, double lowM=100, double highM=1000){

  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // get workspace

  char fileName[100];  
  sprintf(fileName,"/scratch0/hep/whitbeck/4lHelicity/Combination/2012higgsReview2D/workspaces/hzz4l_%sS.%i.0.input.root",channel,mass);
  
  cout << fileName << endl;

  TFile* wspFile = new TFile(fileName);  
  RooWorkspace* w = (RooWorkspace*) wspFile->Get("w");

  if(!w) 
    return;

  w->var("CMS_zz4l_mass")->SetName("zzmass");

  cout << "changed names" << w->var("zzmass") << endl;

  RooRealVar MC_weight("MC_weight","MC_weight",0,1000);

  RooArgSet obs(*(w->var("zzmass")),*(w->var("melaLD")));
  RooArgSet obsW(*(w->var("zzmass")),*(w->var("melaLD")),MC_weight);
  
  cout << " initialized argset " << endl;

  // load MC 
  if(strcmp(channel,"2e2mu")==0)
    channel="2mu2e";

  sprintf(fileName,"../datafiles/ZZ%sAnalysisTree_H%i_withDiscriminants.root",channel,mass);
  cout << fileName << endl;
  TChain* treeMC = new TChain("angles");
 
  treeMC->Add(fileName);
  
   if(!treeMC || treeMC->GetEntries()<=0 )
    return ;

   cout << treeMC->GetEntries() << endl;

   char cutString[100]="";
   sprintf(cutString,"zzmass>%i&&zzmass<%i",(int)lowM,(int)highM);

   RooDataSet* dataMC = new RooDataSet("dataMC","dataMC",treeMC,obsW,cutString,"MC_weight");  // need to first load workspace and change names
   cout << "dataset size: " << dataMC->sumEntries() << endl;

   // make 2D model
   
   TFile* tempFile = new TFile(tempFileName);
   TH2F* temp = (TH2F*) tempFile->Get("h_mzzD");
   TH2F* oldTemp = (TH2F*) tempFile->Get("oldTemp");
   
   cout << "loaded templates " << temp << " " << oldTemp << endl;
   
   RooDataHist* sigTempDataHist = new RooDataHist("sigTempDataHist","sigTempDataHist",obs,temp);
   RooDataHist* oldSigTempDataHist = new RooDataHist("oldSigTempDataHist","oldSigTempDataHist",obs,oldTemp);
  
   cout << "initialized RooDataHist" << endl;
   
   RooHistPdf* sigTemplatePdf_ggH = new RooHistPdf("sigTemplatePdf_ggH","sigTemplatePdf_ggH",obs,*sigTempDataHist);
   RooHistPdf* oldSigTemplatePdf_ggH = new RooHistPdf("oldSigTemplatePdf_ggH","oldSigTemplatePdf_ggH",obs,*sigTempDataHist);
   
   cout << "initialized RooHistPdf " << sigTemplatePdf_ggH << " " << oldSigTemplatePdf_ggH << endl;
   
   RooProdPdf *sig2d_ggH;
   RooProdPdf *oldSig2d_ggH;
   if(mass<170){
     sig2d_ggH =  new RooProdPdf("sig2d_ggH","sig2d_ggH",*(w->pdf("signalCB_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
     oldSig2d_ggH =  new RooProdPdf("oldSig2d_ggH","oldSig2d_ggH",*(w->pdf("signalCB_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
   }
   else{
     sig2d_ggH =  new RooProdPdf("sig2d_ggH","sig2d_ggH",*(w->pdf("sig_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
     oldSig2d_ggH =  new RooProdPdf("oldSig2d_ggH","oldSig2d_ggH",*(w->pdf("sig_ggH")),RooFit::Conditional(*sigTemplatePdf_ggH,RooArgSet(*w->var("melaLD"))));
   }
   
   cout << "initialized 2D PDFs" << endl;
   
   RooDataSet *sigPDF = sig2d_ggH->generate(obs,100000);
   sigPDF->reduce(cutString);
   RooDataSet *oldSigPDF = oldSig2d_ggH->generate(obs,100000);
   sigPDF->reduce(cutString);
   
   // plot MELA distribution
   
   RooPlot* melaPlot = w->var("melaLD")->frame(30);
   dataMC->plotOn(melaPlot,Rescale(sigPDF->sumEntries()/dataMC->sumEntries()));
   sigPDF->plotOn(melaPlot,MarkerStyle(2),MarkerColor(4));
   oldSigPDF->plotOn(melaPlot,MarkerStyle(3),MarkerColor(2));
   //sig2d_ggH->plotOn(melaPlot);
   
   TCanvas* can = new TCanvas("can","can",400,400);
   melaPlot->Draw();
   sprintf(fileName,"2DpdfProjxCheck_mH%i_%s_%i-%i.png",mass,channel,(int)lowM,(int)highM);
   can->SaveAs(fileName);
   sprintf(fileName,"2DpdfProjxCheck_mH%i_%s_%i-%i.eps",mass,channel,(int)lowM,(int)highM);
   can->SaveAs(fileName);
}

void runAllcocktailSyst(char* channel="4mu"){

  int mH[15]={120,130,140,150,160,170,180,190,200,250,300,350,400,500,600};
  
  char temp[100];
  
  sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  
  for(int i=0; i<15; i++){
      cocktailSyst2(mH[i],channel,temp);    
  }

}

// = = = = = = = = = = = = = = = = = = = = 
// depricated
// = = = = = = = = = = = = = = = = = = = = 

void cocktailSyst(int index,char* channel){


  int mH[15]={120,130,140,150,160,170,180,190,200,250,300,350,400,500,600};
  double wH[15]={1.,1.,1.,1.,1.,1.,1.,1.,1.43,4.04,8.43,15.2,29.2,68.,123.};

  char temp[150];
  sprintf(temp,"../datafiles/ZZ%sAnalysisTree_H%i_withDiscriminants.root",channel,mH[index]);

  TChain* tree = new TChain("angles");
  tree->Add(temp);

  if(!tree || tree->GetEntries()<=0 )
    return;

  sprintf(temp,"MC_weight*(zzmass>%i&&zzmass<%i)",int(mH[index]-wH[index]*2),int(mH[index]+wH[index]*2));
  
  cout << "drawing tree: " << temp << endl;
  tree->Draw("melaLD>>hMC(30,0,1)",temp);
  TH1F* histoMC = (TH1F*) gDirectory->Get("hMC");
  histoMC->SetLineColor(4);
  histoMC->SetLineStyle(2);
  histoMC->SetLineWidth(2);

  sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
  cout << "loading template: " << temp << endl;
  TFile* tempFile = new TFile(temp);
  TH2F* Htemp = (TH2F*) tempFile->Get("h_mzzD");

  cout << "lower index: " << int(mH[index]/2.-wH[index]-50+1) << " higher index: " << int(mH[index]/2.+wH[index]-50+1)<< endl;

  TH1F* hTemp = (TH1F*) Htemp->ProjectionY("hTemp",int(mH[index]/2.-wH[index]-50+1),int(mH[index]/2.+wH[index]-50+1));
  hTemp->SetLineColor(2);
  hTemp->SetLineWidth(2);
  hTemp->GetXaxis()->SetTitle("D");
  hTemp->GetYaxis()->SetRangeUser(0,hTemp->GetMaximum()>histoMC->GetMaximum()?1.3*hTemp->GetMaximum():histoMC->GetMaximum()*1.3);

  TCanvas* can = new TCanvas("can","can",400,400);

  hTemp->DrawNormalized();
  histoMC->DrawNormalized("SAME");

  TLegend* leg = new TLegend(.2,.7,.5,.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  
  leg->AddEntry(hTemp,"cocktail","l");
  sprintf(temp,"H%i",mH[index]);
  leg->AddEntry(histoMC,temp,"l");
  leg->Draw();

  sprintf(temp,"cocktailSyst_mH%i_%s.eps",mH[index],channel);
  can->SaveAs(temp);
  sprintf(temp,"cocktailSyst_mH%i_%s.png",mH[index],channel);
  can->SaveAs(temp);
}

void crossCheckSmoothSlices(char* channel="4mu", bool isSig=true, int start=340, int end=350){

  char fileName[100];
  if(isSig)
    sprintf(fileName,"../datafiles/Dsignal_%s.root",channel);
  else
    sprintf(fileName,"../datafiles/Dbackground_%s.root",channel);

  TFile* file = new TFile(fileName);

  TH2F* oldTemp = (TH2F*) file->Get("oldTemp");
  TH2F* newTemp = (TH2F*) file->Get("h_mzzD");
  
  TH1F* oldHisto = (TH1F*) oldTemp->ProjectionY("old",start,end);
  oldHisto->SetLineColor(4);
  oldHisto->SetLineWidth(2);
  TH1F* newHisto = (TH1F*) newTemp->ProjectionY("new",start,end);
  newHisto->SetLineColor(2);
  newHisto->SetLineWidth(2);
  newHisto->SetLineStyle(2);

  oldHisto->Draw();
  newHisto->Draw("SAME");

}

void compareZplusX_SSdata_vs_OSdata(double mzzLow=130, double mzzHigh=180){

  // = = = = = = = = 
  // load trees 
  // = = = = = = = =
  
  TChain* CRssData = new TChain("angles");
  CRssData->Add("../datafiles/ssCR/ZZ*AnalysisTree_Double*_4670_withDiscriminants.root");

  TChain* CRosData = new TChain("angles");
  CRosData->Add("../datafiles/osCR/ZZ*AnalysisTree_Double*_4670_withDiscriminants.root");

  if( !CRosData || CRosData->GetEntries()<=0 || !CRssData || CRssData->GetEntries()<=0 ){
    cout << "problem loading files... " << endl;
    return;
  }

  // = = = = = = =
  // set branches 
  // = = = = = = =

  double mzz, D;

  CRosData->SetBranchAddress("zzmass",&mzz);
  CRosData->SetBranchAddress("melaLD",&D);

  CRssData->SetBranchAddress("zzmass",&mzz);
  CRssData->SetBranchAddress("melaLD",&D);

  // = = = = = = = = = = = 
  // initialize histograms 
  // = = = = = = = = = = = 

  TH1F* h_CRosData = new TH1F("h_CRosData",";MELA;",30,0,1);
  h_CRosData->Sumw2();
  TH1F* h_CRssData = new TH1F("h_CRssData",";MELA;",30,0,1);
  h_CRssData->Sumw2();

  TCanvas* can = new TCanvas("can","can",400,400);

  // = = = = = =
  // fill histos
  // = = = = = = 

  for(int iEvt=0; iEvt<CRosData->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRosData->GetEntries() <<endl;
    CRosData->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRosData->Fill(D);
    }
    
  }

  for(int iEvt=0; iEvt<CRssData->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRssData->GetEntries() <<endl;
    CRssData->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRssData->Fill(D);
    }
    
  }
  
  // = = = = = = 
  // draw histos
  // = = = = = = 

  h_CRosData->Scale(h_CRssData->Integral()/h_CRosData->Integral());
  h_CRosData->SetLineColor(2);
  h_CRosData->SetLineWidth(2);

  h_CRssData->SetMarkerStyle(8);

  if(h_CRssData->GetMaximum()+h_CRssData->GetBinError(h_CRssData->GetMaximumBin()) > h_CRosData->GetMaximum()+h_CRosData->GetBinError(h_CRosData->GetMaximumBin())){
    h_CRssData->GetYaxis()->SetRangeUser(0,(h_CRssData->GetMaximum()+h_CRssData->GetBinError(h_CRssData->GetMaximumBin()))*1.5);
    h_CRssData->Draw("ep");
    h_CRosData->Draw("EHISTSAME");
  }else{
    h_CRosData->GetYaxis()->SetRangeUser(0,(h_CRosData->GetMaximum()+h_CRosData->GetBinError(h_CRosData->GetMaximumBin()))*1.5);
    h_CRosData->Draw("EHIST");
    h_CRssData->Draw("SAMEep");
  }

  // ---------- LEGEND ---------------

  TLegend* leg = new TLegend(.5,.6,.90,.90);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_CRssData,"Z+X 2P+2F SS (data)","lp");
  leg->AddEntry(h_CRosData,"Z+X 2P+2F OS (data)","l");
  //leg->AddEntry(h_qqZZ,"qqZZ (Powheg MC)","l");
  
  leg->Draw();


  char temp[50];
  sprintf(temp,"Z+X_SSdata_vs_OSdata_%i-%i.eps",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);
  sprintf(temp,"Z+X_SSdata_vs_OSdata_%i-%i.png",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);

}

void runSSvsOS(){
  
  compareZplusX_SSdata_vs_OSdata(100,110);
  compareZplusX_SSdata_vs_OSdata(110,120);
  compareZplusX_SSdata_vs_OSdata(120,130);
  compareZplusX_SSdata_vs_OSdata(130,140);
  compareZplusX_SSdata_vs_OSdata(140,150);
  compareZplusX_SSdata_vs_OSdata(150,160);
  compareZplusX_SSdata_vs_OSdata(160,170);
  compareZplusX_SSdata_vs_OSdata(170,180);
    
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// compare heavy flavor MC with light flavor mC  in SS+OS CR
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

void compareZplusX_Bs_vs_noBs(double mzzLow=130, double mzzHigh=180){

  // = = = = = = = = 
  // load trees 
  // = = = = = = = =
  
  TChain* CRwithBs = new TChain("angles");
  CRwithBs->Add("../datafiles/summedCR/ZZ*AnalysisTree_DYJetsToLLTuneZ2*_B_withDiscriminants.root");

  TChain* CRnoBs = new TChain("angles");
  CRnoBs->Add("../datafiles/summedCR/ZZ*AnalysisTree_DYJetsToLLTuneZ2*_NoB_withDiscriminants.root");

  if( !CRnoBs || CRnoBs->GetEntries()<=0 || !CRwithBs || CRwithBs->GetEntries()<=0 ){
    cout << "problem loading files... " << endl;
    return;
  }

  // -------------
  // set branches 
  // -------------

  double mzz, D;

  CRnoBs->SetBranchAddress("zzmass",&mzz);
  CRnoBs->SetBranchAddress("melaLD",&D);

  CRwithBs->SetBranchAddress("zzmass",&mzz);
  CRwithBs->SetBranchAddress("melaLD",&D);

  // ---------------------
  // initialize histograms 
  // ---------------------

  TH1F* h_CRnoBs = new TH1F("h_CRnoBs",";MELA;",30,0,1);
  h_CRnoBs->Sumw2();
  TH1F* h_CRwithBs = new TH1F("h_CRwithBs",";MELA;",30,0,1);
  h_CRwithBs->Sumw2();

  TCanvas* can = new TCanvas("can","can",400,400);

  // ------------
  // fill histos
  // ------------

  for(int iEvt=0; iEvt<CRnoBs->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRnoBs->GetEntries() <<endl;
    CRnoBs->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRnoBs->Fill(D);
    }
    
  }

  for(int iEvt=0; iEvt<CRwithBs->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRwithBs->GetEntries() <<endl;
    CRwithBs->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRwithBs->Fill(D);
    }
    
  }
  
  // -----------
  // draw histos
  // -----------

  h_CRnoBs->Scale(h_CRwithBs->Integral()/h_CRnoBs->Integral());
  h_CRnoBs->SetLineColor(2);
  h_CRnoBs->SetLineWidth(2);

  h_CRwithBs->SetMarkerStyle(8);

  if(h_CRwithBs->GetMaximum()+h_CRwithBs->GetBinError(h_CRwithBs->GetMaximumBin()) > h_CRnoBs->GetMaximum()+h_CRnoBs->GetBinError(h_CRnoBs->GetMaximumBin())){
    h_CRwithBs->GetYaxis()->SetRangeUser(0,(h_CRwithBs->GetMaximum()+h_CRwithBs->GetBinError(h_CRwithBs->GetMaximumBin()))*1.5);
    h_CRwithBs->Draw("ep");
    h_CRnoBs->Draw("EHISTSAME");
  }else{
    h_CRnoBs->GetYaxis()->SetRangeUser(0,(h_CRnoBs->GetMaximum()+h_CRnoBs->GetBinError(h_CRnoBs->GetMaximumBin()))*1.5);
    h_CRnoBs->Draw("EHIST");
    h_CRwithBs->Draw("SAMEep");
  }

  // ---------- LEGEND ---------------

  TLegend* leg = new TLegend(.5,.6,.90,.90);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_CRwithBs,"Z+X 2P+2F HF (MadGraph MC)","lp");
  leg->AddEntry(h_CRnoBs,"Z+X 2P+2F LF (MadGraph MC)","l");
  //leg->AddEntry(h_qqZZ,"qqZZ (Powheg MC)","l");
  
  leg->Draw();


  char temp[50];
  sprintf(temp,"Z+X_MCwithBs_vs_MCnoBs_%i-%i.eps",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);
  sprintf(temp,"Z+X_MCwithBs_vs_MCnoBs_%i-%i.png",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);

}

void runBvsNoBs(){

  compareZplusX_Bs_vs_noBs(100,110);
  compareZplusX_Bs_vs_noBs(110,120);
  compareZplusX_Bs_vs_noBs(120,130);
  compareZplusX_Bs_vs_noBs(130,140);
  compareZplusX_Bs_vs_noBs(140,150);
  compareZplusX_Bs_vs_noBs(150,160);
  compareZplusX_Bs_vs_noBs(160,170);
  compareZplusX_Bs_vs_noBs(170,180);
    
}
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// function below plots qqZZ, data CR, MC CR shapes, calculates
// the ratio and fits the ratio with a line
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

void measureSystematics(double mzzLow=130, double mzzHigh=180){

  // ----------
  // load trees 
  // ----------
  
  TChain* CRdata = new TChain("SelectedTree");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CREEEEssTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CRMMMMssTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CREEMMssTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CRMMEEssTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CREEEEosTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CRMMMMosTree.root");
  CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CREEMMosTree.root");
  //CRdata->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_Double*_CRMMEEosTree.root");

  TChain* CRmc = new TChain("SelectedTree");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CREEEEssTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRMMMMssTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CREEMMssTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRMMEEssTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CREEEEosTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRMMMMosTree.root");
  CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CREEMMosTree.root");
  //CRmc->Add("CJLSTtrees_June21_2012/7TeV_FSR/CR/HZZ4lTree_DYJetsToLLTuneZ2*_CRMMEEosTree.root");

  TChain* qqZZ = new TChain("SelectedTree");
  qqZZ->Add("CJLSTtrees_June21_2012/7TeV_FSR/HZZ*Tree_ZZTo*.root");

  if( !CRmc || CRmc->GetEntries()<=0 || !CRdata || CRdata->GetEntries()<=0 || !qqZZ || qqZZ->GetEntries()<=0 ){
    cout << "problem loading files... " << endl;
    return;
  }

  // = = = = = = =
  // set branches 
  // = = = = = = =

  float mzz, D, w;

  CRmc->SetBranchAddress("ZZMass",&mzz);
  CRmc->SetBranchAddress("ZZLD",&D);
  CRmc->SetBranchAddress("MC_weight",&w);

  CRdata->SetBranchAddress("ZZMass",&mzz);
  CRdata->SetBranchAddress("ZZLD",&D);

  
  qqZZ->SetBranchAddress("ZZMass",&mzz);
  qqZZ->SetBranchAddress("ZZLD",&D);
  qqZZ->SetBranchAddress("MC_weight",&w);

  // = = = = = = = = = = = 
  // initialize histograms 
  // = = = = = = = = = = = 

  TH1F* h_CRmc = new TH1F("h_CRmc",";MELA;",30,0,1);
  h_CRmc->Sumw2();
  TH1F* h_CRdata = new TH1F("h_CRdata",";MELA;",30,0,1);
  h_CRdata->Sumw2();
  TH1F* h_qqZZ = new TH1F("h_qqZZ",";MELA;",30,0,1);
  h_qqZZ->Sumw2();

  TCanvas* can = new TCanvas("can","can",400,550);
  TPad* pad2 = new TPad("pad2","pad2",0.,0.,1.,0.3);
  TPad* pad1 = new TPad("pad1","pad1",0.,0.3,1.,1.);
  pad2->SetBottomMargin(0.22);
  
  pad1->Draw();
  pad2->Draw();

  // = = = = = =
  // fill histos
  // = = = = = = 


  for(int iEvt=0; iEvt<qqZZ->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << qqZZ->GetEntries() <<endl;
    qqZZ->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_qqZZ->Fill(D,w);
    }
    
  }

  for(int iEvt=0; iEvt<CRmc->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRmc->GetEntries() <<endl;
    CRmc->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRmc->Fill(D,w);
    }
    
  }

  for(int iEvt=0; iEvt<CRdata->GetEntries(); iEvt++){
    
    //if(iEvt%1000 == 0) cout << "event " << iEvt << "/" << CRdata->GetEntries() <<endl;
    CRdata->GetEntry(iEvt);

    if(mzz>mzzLow&&mzz<mzzHigh){
      h_CRdata->Fill(D);
    }
    
  }
  
  cout << "about to draw" << endl;

  // = = = = = = 
  // draw histos
  // = = = = = = 

  h_qqZZ->Scale(h_CRdata->Integral()/h_qqZZ->Integral());
  h_qqZZ->SetLineColor(4);
  h_qqZZ->SetLineStyle(2);  
  h_qqZZ->SetLineWidth(2);

  cout << "check" <<endl;

  h_CRmc->Scale(h_CRdata->Integral()/h_CRmc->Integral());
  h_CRmc->SetLineColor(2);
  h_CRmc->SetLineWidth(2);

  h_CRdata->SetMarkerStyle(8);

  cout << "style set!" << endl;

  pad1->cd();

  if(h_CRdata->GetMaximum()+h_CRdata->GetBinError(h_CRdata->GetMaximumBin()) > h_CRmc->GetMaximum()+h_CRmc->GetBinError(h_CRmc->GetMaximumBin())){
    h_CRdata->GetYaxis()->SetRangeUser(0,(h_CRdata->GetMaximum()+h_CRdata->GetBinError(h_CRdata->GetMaximumBin()))*1.5);
    h_CRdata->Draw("ep");
    h_CRmc->Draw("EhistSAME");
    h_qqZZ->Draw("EhistSAME");
  }else{
    h_CRmc->GetYaxis()->SetRangeUser(0,(h_CRmc->GetMaximum()+h_CRmc->GetBinError(h_CRmc->GetMaximumBin()))*1.5);
    h_CRmc->Draw("Ehist");
    h_CRdata->Draw("SAMEep");
    h_qqZZ->Draw("EhistSAME");
  }

  cout << "done drawing?" << endl;

  // ---------- LEGEND ---------------

  TLegend* leg = new TLegend(.5,.6,.90,.90);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_CRdata,"Z+X 2P+2F (data)","lp");
  leg->AddEntry(h_CRmc,"Z+X 2P+2F (Madgraph MC)","l");
  leg->AddEntry(h_qqZZ,"qqZZ (Powheg MC)","l");
  
  leg->Draw();

  // ------------ RATIO PAD ---------------

  cout << "ratio" << endl;

  pad2->cd();
  
  TH1F* ratio_qqZZ   = new TH1F(*h_qqZZ);
  ratio_qqZZ->Divide(h_qqZZ);
  TH1F* ratio_CRmc   = new TH1F(*h_CRmc);
  ratio_CRmc->Divide(h_qqZZ);
  TH1F* ratio_CRdata = new TH1F(*h_CRdata);
  ratio_CRdata->Divide(h_qqZZ);

  ratio_CRdata->GetYaxis()->SetRangeUser(0.,2.);
  ratio_CRdata->Draw("p");

  // ------------- FIT RATIO -------------

  cout << "Fitting data CR " << mzzLow << "-" << mzzHigh << endl;

  TF1* fline = new TF1("fline","[0]+[1]*x",mzzLow,mzzHigh);
  fline->SetLineWidth(2);
  fline->SetLineColor(kGreen+1);
  ratio_CRdata->Fit("fline","");
  //gStyle->SetOptFit(0);

  cout << " done with fits drawing CRmc and qqZZ " << endl;

  ratio_CRmc->Draw("EhistSAME");
  ratio_qqZZ->Draw("EhistSAME");

  // -------------------------------------

  char temp[50];
  sprintf(temp,"Z+X_data_vs_MC_vs_qqZZ_%i-%i.eps",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);
  sprintf(temp,"Z+X_data_vs_MC_vs_qqZZ_%i-%i.png",(int)mzzLow,(int)mzzHigh);
  can->SaveAs(temp);

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// measure systematic effect in course binned mZZ windows
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

void runSystematics(){

  measureSystematics(100,110);
  measureSystematics(110,120);
  measureSystematics(120,130);
  measureSystematics(130,140);
  measureSystematics(140,150);
  measureSystematics(150,160);
  measureSystematics(160,170);
  measureSystematics(170,180);
  measureSystematics(100,180);

  /*
  measureSystematics(100,120);
  measureSystematics(120,140);
  measureSystematics(140,160);
  measureSystematics(160,180);
  measureSystematics(180,220);
  measureSystematics(220,260);
  measureSystematics(260,300);
  */


}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
// compare channels slice by slice
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

void channelCompareSlice(char* sample="signal"){

  char temp[100];

  sprintf(temp,"../datafiles/D%s_4mu.root",sample);
  TFile* f4mu = new TFile(temp);
  TH2F* tmp4mu = (TH2F*) f4mu->Get("oldTemp");
  TH1F* h4mu;

  sprintf(temp,"../datafiles/D%s_4e.root",sample);
  TFile* f4e = new TFile(temp);
  TH2F* tmp4e = (TH2F*) f4e->Get("oldTemp");
  TH1F* h4e;

  sprintf(temp,"../datafiles/D%s_2e2mu.root",sample);
  TFile* f2e2mu = new TFile(temp);
  TH2F* tmp2e2mu = (TH2F*) f2e2mu->Get("oldTemp");
  TH1F* h2e2mu;

  sprintf(temp,"channelComparison_%s",sample);
  TCanvas* can = new TCanvas(temp,temp,400,400);


  for(int i=41;i<108;i++){

    
    h4mu = (TH1F*) tmp4mu->ProjectionY("4mu",i,i);
    h4mu->SetLineColor(1);
    h4mu->Draw("");

    h4e = (TH1F*) tmp4e->ProjectionY("4e",i,i);
    h4e->SetLineColor(2);
    h4e->Draw("SAME");

    h2e2mu = (TH1F*) tmp2e2mu->ProjectionY("2e2mu",i,i);
    h2e2mu->SetLineColor(4);
    h2e2mu->Draw("SAME");
   
    if(i==41)
      can->Print(".pdf(","pdf");
    else 
      can->Print(".pdf","pdf");
  }
  can->Print(".pdf)","pdf");


}

