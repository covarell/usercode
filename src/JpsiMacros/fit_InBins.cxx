 // C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>

#include <TLatex.h>
#include "RooFitResult.h"
// ROOT includes
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2F.h>
#include <TH1F.h>
#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
//#include "CommonVar.h"
using namespace RooFit;



void defineBackground(RooWorkspace *ws)
{
  //BACKGROUND
  cout << "defining background" << endl;
  //ws->factory("Polynomial::CPolFunct(Jpsi_Mass,{CoefPol1[-1.0,-300.,1000.]})");
  //ws->factory("Polynomial::CPolFunct2(Jpsi_Mass,{CoefPol1,CoefPol2[0.5,0.,1000]})");
  //ws->factory("Chebychev::expFunct(Jpsi_Mass,{coefExp[0.,-1.,1],coefExpB[0.,0.,0.]})");
  ws->factory("Exponential::expFunct(Jpsi_Mass,coefExp[-5,-10.,4])");

  return;
}

void defineSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:
  cout << "defining signal" << endl;
  //Normal Gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.1,3.05,3.15],sigmaSig1[0.03,0.008,0.08])");
  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig1,sigmaSig2[0.04,0.008,0.15])");
  // ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig1,sigmaSig1");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(Jpsi_Mass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(Jpsi_Mass,meanSig1,sigmaSig1,alpha[0.5,.1,3.],enne[8.,2.,10.])");

 ws->factory("CBShape::sigCB2(Jpsi_Mass,meanSig1,sigmaSigCB2[0.04,0.,4.],alpha1[0.5,0.,4.],enne1[5.,1.,50.])");
  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.01,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  //ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");
  ws->factory("SUM::sigCBGauss(coeffGauss*signalG1,signalG2)");

  //Sum of Two CBs
  ws->factory("SUM::sigDoubleCB(coeffGauss*sigCB,sigCB2)");

  return;
}

void prefitSideband(RooWorkspace *ws, const int DataCat, int k, int j)
{
  cout << "prefit Sideband" << endl;
  ws->var("coefExp")->setConstant(kFALSE);
  ws->var("alpha")->setConstant(kFALSE);
  ws->var("enne")->setConstant(kFALSE);
  // ws->var("CoefPol1")->setConstant(kFALSE);
  //ws->var("CoefPol2")->setConstant(kFALSE);
  ws->var("sigmaSig2")->setConstant(kFALSE);
  ws->var("meanSig1")->setConstant(kFALSE);
  RooDataSet *tmpdata;
  tmpdata = (RooDataSet*)ws->data("data");
  ws->pdf("expFunct")->fitTo(*tmpdata,Range("left,right"),SumW2Error(kTRUE));
  //ws->pdf("sigCBGauss")->fitTo(*tmpdata,Range("peak"),SumW2Error(kTRUE));
  // if (k==3 && (j==13 || j==7 || j==3 || j==14)){
  ws->var("alpha")->setVal(1.8);
  ws->var("alpha")->setConstant(kTRUE);
  //ws->var("enne")->setVal(4);
  //ws->var("enne")->setConstant(kTRUE);
  //}
  // ws->var("coefExp")->setConstant(kTRUE);
  //ws->var("alpha")->setVal(1.5);
  //ws->var("alpha")->setConstant(kTRUE);
  //ws->var("enne")->setVal(4);
  //ws->var("enne")->setConstant(kTRUE);
  // ws->var("CoefPol1")->setConstant(kTRUE);
  // ws->var("CoefPol2")->setConstant(kTRUE);
  // ws->var("sigmaSig2")->setConstant(kTRUE);
  //ws->var("meanSig1")->setConstant(kTRUE);
  
  return;
}

void setRanges(RooWorkspace *ws)
{
  cout<< "setting ranges" << endl;
  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.4;

  ws->var("Jpsi_Mass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("Jpsi_Mass")->setRange("left1",2.75,2.9);
  ws->var("Jpsi_Mass")->setRange("left",JpsiMassMin,2.95);
  ws->var("Jpsi_Mass")->setRange("right",3.25,JpsiMassMax);
  ws->var("Jpsi_Mass")->setRange("right1",3.4,JpsiMassMax);
  ws->var("Jpsi_Mass")->setRange("peak",2.95,3.15);

  ws->cat("Jpsi_Type")->setRange("glbglb","GG");
  ws->cat("Jpsi_Type")->setRange("glbtrk","GT");
  ws->cat("Jpsi_Type")->setRange("trktrk","TT");

  return;
}

double drawResults(RooWorkspace *ws, const int npar, char *hst, TH2F * hrms, int bin, double pt1, double pt2, double eta1, double eta2)

{
  cout << "drawing results" << endl;
 
  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");

  char reducestr[200];
  const double sigma = ws->var("sigmaSig1")->getVal();
  const double mean = ws->var("meanSig1")->getVal();
  double up=mean+(2.5*sigma);
  double down=mean-(2.5*sigma);

  RooPlot *mframe = ws->var("Jpsi_Mass")->frame(2.6,3.4);
  RooPlot *ptframe = ws->var("Jpsi_Pt")->frame();
  RooPlot *etaframe = ws->var("Jpsi_Y")->frame();


  mframe->SetTitle("");
  
  //Plotting Data    
  data->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(40),MarkerStyle(7),MarkerSize(.5));
  data->statOn(ptframe);
  data->plotOn(ptframe,DataError(RooAbsData::SumW2),Cut("Jpsi_Mass>= 3.0 && Jpsi_Mass <= 3.2"),Binning(100));
  data->plotOn(etaframe,DataError(RooAbsData::SumW2),Cut("Jpsi_Mass>= 3.0 && Jpsi_Mass <= 3.2"),Binning(100));
  RooDataSet* datar=(RooDataSet*)data->reduce("Jpsi_Mass>= 3.0 && Jpsi_Mass <= 3.2");
  RooRealVar *meanv= datar->meanVar(*ws->var("Jpsi_Pt"));
  RooRealVar *rmsv= datar->rmsVar(*ws->var("Jpsi_Pt"));
  hrms->SetBinContent(bin,meanv->getVal());
  hrms->SetBinError(bin,rmsv->getVal());
  
  double pmax=0, pmin=100;
  if (data->sumEntries()!=0){
    TH1F* hd= new TH1F("hd","hd",40,2.6,3.4); 
    data->fillHistogram(hd,RooArgList(RooArgSet(*ws->var("Jpsi_Mass"))));
    //hd= dynamic_cast<TH1F*>(data->createHistogram("histod",*ws->var("Jpsi_Mass"),AutoBinning(160)));
    pmax=hd->GetMaximum();
    pmin=100;
    double pmintmp=10000;
    for (int s=1; s!=hd->GetNbinsX(); s++){
      if (hd->GetBinContent(s)!=0){
	if (hd->GetBinContent(s) < pmintmp){
	  pmin=hd->GetBinContent(s);
	  pmintmp=pmin;
	}
      }
    }
    delete hd;
  }
  
  //Plotting PDF
  totPDF->plotOn(mframe,LineColor(kBlack),LineWidth(1.0),Normalization(1.0,RooAbsReal::RelativeExpected));
  totPDF->paramOn(mframe,Layout(0.7,0.95,0.95),Parameters(RooArgSet(*ws->var("NSig"),*ws->var("NBkg"))));
  double chi=mframe->chiSquare(npar);
  
  totPDF->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineWidth(1.0),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));
  
  mframe->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");

  if (pmax>pmin) {
    mframe->GetYaxis()->SetRangeUser(0.6*pmin,pmax*3); 
    if (pmin==0) mframe->GetYaxis()->SetRangeUser(0.5,pmax*3);
  }

  TCanvas c1;
  
  c1.cd();
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(13);
  t->SetTextFont(63);
  t->SetTextSizePixels(25);
  c1.SetLogy(1);
  char chisq[40];
  mframe->Draw();
  
  sprintf(chisq,"%5.2f #leq P_{T} < %5.2f",pt1,pt2);
  t->DrawLatex(0.2,0.90,chisq);
  sprintf(chisq,"%5.2f #leq Y < %5.2f",eta1,eta2);
  t->DrawLatex(0.2,0.85,chisq);
  //sprintf(chisq,"#chi^{2} = %4.2f",chi);
  //t->DrawLatex(0.2,0.80,chisq);
  
  //c1.Modified();
  //c1.Update();
  
  
  sprintf(reducestr,"InBins/massfit_%s.pdf",hst);
  c1.SaveAs(reducestr);
  sprintf(reducestr,"InBins/massfit_%s.png",hst);
  c1.SaveAs(reducestr);
  sprintf(reducestr,"InBins/massfit_%s.root",hst);
  c1.SaveAs(reducestr);
  
  sprintf(reducestr,"InBins/ptplot_%s.png",hst);
  TCanvas c2;
  c2.cd();
  ptframe->Draw();
  c2.SaveAs(reducestr);
  
  sprintf(reducestr,"InBins/etaplot_%s.png",hst);
  TCanvas c3;
  c3.cd();
  etaframe->Draw();
  c3.SaveAs(reducestr);
  
  delete t;
  return chi;
}

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol,double &Mean, double &errMean, double &ErrResol, double &SN, double &Bkg, TH2F * h,int bin, RooFitResult * rfr,int npar,char* hst, double* pt_bins, double* y_bins, int SIZEY, int SIZEPT, double JPT, double JY)
{
  cout << "printing " << endl;
  Nsig   = ws->var("NSig")->getVal();
  const double Nbkg = ws->var("NBkg")->getVal();
  errSig = ws->var("NSig")->getError();
  const double coeffGauss = ws->var("coeffGauss")->getVal();
  const double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  //ErrResol=ws->var("sigmaSig1")->getError();
  Mean = ws->var("meanSig1")->getVal();
  errMean=ws->var("meanSig1")->getError();
  const double sigmaSig2 = ws->var("sigmaSig2")->getVal();
  const double ecoeffGauss = ws->var("coeffGauss")->getError();
  const double esigmaSig1 = ws->var("sigmaSig1")->getError();
  const double esigmaSig2 = ws->var("sigmaSig2")->getError();
  resol = sigmaSig1;
  ErrResol = esigmaSig1;
  //resol = sqrt(coeffGauss*sigmaSig1*sigmaSig1 + (1-coeffGauss)*sigmaSig2*sigmaSig2);
  //ErrResol = (0.5/resol)*sqrt(pow(sigmaSig1*coeffGauss*esigmaSig1,2) + pow(sigmaSig2*(1-coeffGauss)*esigmaSig2,2) + pow(0.5*(sigmaSig1*sigmaSig1 - sigmaSig2*sigmaSig2)*ecoeffGauss,2));
  
  ws->var("Jpsi_Mass")->setRange("integral",Mean-2.5*resol,Mean+2.5*resol);
  
  RooAbsReal* N= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("integral"));
  //RooAbsReal* S= ws->pdf("sigCBGauss")->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("integral"));
  RooAbsReal* S= ws->pdf("sigCB")->createIntegral(RooArgSet(*ws->var("Jpsi_Mass")),NormSet(RooArgSet(*ws->var("Jpsi_Mass"))),Range("integral"));  
  
  Bkg=Nbkg*N->getVal();
  
  SN=(Nsig*S->getVal())/(Nbkg*N->getVal());
  h->SetBinContent(bin,Nsig);
  h->SetBinError(bin,errSig);
 
  char hname[50];

  TH2F * HresultsJpsiBin= new TH2F("YieldsJBin","YieldsJBin",SIZEY,y_bins,SIZEPT,pt_bins);
 
  HresultsJpsiBin->SetBinContent(bin,Nsig);
  HresultsJpsiBin->SetBinError(bin,errSig);

  sprintf(hname,"InBins/HBin_Results_%s.root",hst);

  TFile Binout(hname,"RECREATE");
  Binout.cd();
  HresultsJpsiBin->Write();
  Binout.Close();
  delete HresultsJpsiBin;

  sprintf(hname,"InBins/Results_Bin_%s.txt",hst);
  std::ofstream resultsBin;
  resultsBin.open(hname);
  resultsBin << "-------------------------NEW BIN-----------------------" << endl;
  resultsBin << "Pt bin center " << JPT << "   Y bin center " << JY << endl;
  resultsBin << "Name in plot  " << hst << endl;
  resultsBin << "-------------------------------------------------------" <<endl;
  resultsBin << "Fit : " << Nsig << " +/- " << errSig << endl;
  resultsBin << "Resolution : " << resol*1000. << " MeV  +/- "<<ErrResol*1000  << endl; 
  resultsBin << "Mean : " << Mean << " GeV  +/- "<<errMean  << endl;
  resultsBin << "S/B= " << SN << endl; 
  resultsBin << "-------------------------------------------------------" <<endl;
  resultsBin << "----------------------------------J/psi-------------------------------" << endl;
  rfr->printMultiline(resultsBin,npar);
  resultsBin << "-------------------------------------------------------" <<endl;
  resultsBin << endl;
  resultsBin << endl;  
  resultsBin.close();

  return;
}

int main(int argc, char* argv[])
{
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".L setTDRStyle_modified.C");
  gROOT->ProcessLine("setTDRStyle()");

  char inputFileName[150];
  char *filename;
  bool sidebandPrefit = false;
 
  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];
    
    switch(pchar[0]){
      
    case '-':{
      
      switch(pchar[1]){
	
      case 'f':{
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
      
      case 's':{
        sidebandPrefit = true;
        cout << "Sideband pre-fitting activated" << endl;
        break;
      }
      }
    }
    }
  }

  //REFERENCE BINS {2, 2.75, 3.5, 4.5, 5.5, 6.5, 8, 9, 10, 12, 15, 30}
  Double_t jpsi_y_bins[] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
 
  //Y1
  Double_t jpsi_pt_bins0[] = {7, 8,  9,  10,  11,  12, 13.5,  15, 18, 30, 45, 70};
  //  Double_t jpsi_pt_bins0[] = {7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11.05, 11.35, 11.65, 12, 12.40, 12.80,13.20, 13.7, 14.30, 15, 15.80, 16.80, 18, 20, 23, 30};
  
  //Y2
  Double_t jpsi_pt_bins1[] = {6.5,  8,  9, 10,  12,  15,  30, 45};
 // Double_t jpsi_pt_bins1[] = {6.5, 7, 7.35, 7.7, 8, 8.30, 8.6, 9, 9.5, 10, 10.6, 11.2, 12, 12.8, 13.7, 15, 17.5, 30};

  //Y3
  Double_t jpsi_pt_bins2[] = {4, 5, 6.0,  6.5,  7,  7.5,  8, 8.5, 9,   10,  11,  12,  15,  30, 45};
  
  // Double_t jpsi_pt_bins2[] = {4.5, 5, 5.3, 5.5, 5.7, 5.9, 6.05, 6.20, 6.35, 6.5, 6.65, 6.8, 6.95, 7.1, 7.25, 7.4, 7.55, 7.7, 7.85, 8, 8.15, 8.3, 8.45, 8.60, 8.80, 9,9.2, 9.4, 9.7, 10, 10.3, 10.7, 11.1, 11.5, 12, 12.8, 13.6, 15, 17.5, 30};
  //Y4
  Double_t jpsi_pt_bins3[] = {3.5, 4., 4.25, 4.5, 4.75,  5., 5.25,  5.5, 5.75,  6, 6.25, 6.5, 6.75, 7, 7.25, 7.5,  8, 8.5,  9,  10, 11,  12,  15,  30};

 //Double_t jpsi_pt_bins3[] = {1.5, 2, 2.30, 2.55, 2.75, 2.95, 3.10, 3.25, 3.35, 3.5, 3.65, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,5., 5.1, 5.2, 5.3, 5.4, 5.5, 5.6,5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.45, 7.60, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.25, 9.5, 9.75, 10, 10.5,11,11.5, 12, 13, 15, 17., 30};

  //Y5
  Double_t jpsi_pt_bins4[] =  {2.5, 3, 3.5, 4.,  4.5,  5.5, 6.5,  8,  9,  10, 12, 15, 30};


  int sY=(sizeof(jpsi_y_bins)/sizeof(double))-1;
  int sp0=(sizeof(jpsi_pt_bins0)/sizeof(double))-1;
  int sp1=(sizeof(jpsi_pt_bins1)/sizeof(double))-1;
  int sp2=(sizeof(jpsi_pt_bins2)/sizeof(double))-1;
  int sp3=(sizeof(jpsi_pt_bins3)/sizeof(double))-1;
  int sp4=(sizeof(jpsi_pt_bins4)/sizeof(double))-1;
  cout << "Sizes Y1= " << sp0 << "  Y2= " << sp1 << "  Y3= " << sp2 << "  Y4= " <<sp3 << " Y5= " << sp4 << endl;
  cout << "tot bins= " << sp0+sp1+sp2+sp3+sp4 << endl;
  
  TH2F * Hresults[sY]; //new TH2F("Yields","Yields",sY,jpsi_y_bins,sp,jpsi_pt_bins);
  TH2F * HRMS[sY]; //new TH2F("PtMedioRms","PtMedioRms",sY,jpsi_y_bins,sp,jpsi_pt_bins);
  double JY,JPT;
  
  char name[300];
 
  TFile fIn(filename);
  fIn.cd();
  RooDataSet *reddata2 = (RooDataSet*)fIn.Get("dataJpsi");
  dataJpsi->SetName("data");

  // TFile out("InBins/HistoResults.root","RECREATE");

  for (int k=0; k<sY; k++){
    int sp;
    double *jpsi_pt_bins;
    if (k==0) {jpsi_pt_bins=jpsi_pt_bins0;sp=sp0;}
    if (k==1) {jpsi_pt_bins=jpsi_pt_bins1;sp=sp1;}
    if (k==2) {jpsi_pt_bins=jpsi_pt_bins2;sp=sp2;}
    if (k==3) {jpsi_pt_bins=jpsi_pt_bins3;sp=sp3;}
    if (k==4) {jpsi_pt_bins=jpsi_pt_bins4;sp=sp4;}
    

    char nameh[30];

    sprintf(nameh,"InBins/HistoResults_y%i.root",k+1);
    TFile out(nameh,"RECREATE");

    sprintf(nameh,"Yields_y%i",k+1);
    Hresults[k]= new TH2F(nameh,nameh,sY,jpsi_y_bins,sp,jpsi_pt_bins);
    sprintf(nameh,"PtMedioRms_%i",k+1);
    HRMS[k]= new TH2F(nameh,nameh,sY,jpsi_y_bins,sp,jpsi_pt_bins);

    sprintf(name,"InBins/FitResults_%i.txt",k+1); 
    std::ofstream results;
    results.open(name);

    for (int j=0; j<sp; j++){
      if (k!=2000 /*&& (j==13 || j==7 || j==3 || j==14)*/){
	JY=(jpsi_y_bins[k]+jpsi_y_bins[k+1])/2;
	JPT=(jpsi_pt_bins[j]+jpsi_pt_bins[j+1])/2;

      char cutstring[300];
      char hst[30];
  
      sprintf(cutstring,"((Jpsi_Y > -%5.2f && Jpsi_Y <= -%5.2f) || (Jpsi_Y >= %5.2f && Jpsi_Y < %5.2f)) && (Jpsi_Pt >= %5.2f && Jpsi_Pt < %5.2f)",jpsi_y_bins[k+1],jpsi_y_bins[k],jpsi_y_bins[k],jpsi_y_bins[k+1],jpsi_pt_bins[j],jpsi_pt_bins[j+1]); 
      
      sprintf(hst,"PT%i_Y%i",j+1,k+1);
      cout << " Cut on data: " << cutstring << endl;

      cout << "creating ws"<< endl;
      RooWorkspace *ws = new RooWorkspace("ws");
      cout <<" data reducing" << endl; 
      RooDataSet *reddata= (RooDataSet*)reddata2->reduce(cutstring);
 
      cout << "importing data" << endl;
      ws->import(*reddata);
      cout << " setting ranges" << endl;
      setRanges(ws);

      cout << "######################################## FITTING BIN "<< hst <<" #############################################"<< endl;

      //DEFINE SIGNAL AND BACKGROUND
      defineSignal(ws);
      defineBackground(ws);
     
      ws->factory("SUM::totPDF(NSig[5000.,100.,100000.]*sigCB,NBkg[300.,50.,100000.]*expFunct)");
      int bin=Hresults[k]->FindBin(JY,JPT);

      if (sidebandPrefit) prefitSideband(ws,k,k,j);

      RooFitResult *rfr =ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Save(1),Minos(0),NumCPU(6),SumW2Error(kTRUE));
     
      int npar=rfr->floatParsFinal().getSize();
      double NSigA, errSigA,resolA,meanA,errmeanA,errresolA,SNA,chi,BA;
      chi=drawResults(ws,npar,hst,HRMS[k],bin,jpsi_pt_bins[j],jpsi_pt_bins[j+1],jpsi_y_bins[k],jpsi_y_bins[k+1]);
      printResults(ws,NSigA,errSigA,resolA,meanA,errmeanA,errresolA,SNA,BA,Hresults[k],bin,rfr,npar,hst,jpsi_pt_bins,jpsi_y_bins,sY,sp,JPT,JY);

      cout << "-------------------------NEW BIN----------------------" << endl;
      //cout << "bin definition:  "<< cutstring << endl;
      cout << "Pt in " << jpsi_pt_bins[j] << "," << jpsi_pt_bins[j+1] << endl;
      cout << "Y in " << jpsi_y_bins[k] << "," <<  jpsi_y_bins[k+1]  << endl;
      cout << "Pt bin center " << JPT << "   Y bin center " << JY << endl;
      cout << "Name in plot  " << hst << endl;
      cout << "-------------------------------------------------------" <<endl;
      cout << "Fit : " << NSigA << " +/- " << errSigA << endl;
      cout << "Resolution : " << resolA*1000. << " MeV  +/- "<<errresolA*1000  << endl; 
      cout << "Mean : " << meanA << " GeV  +/- "<<errmeanA  << endl;
      cout << "S/B= " << SNA << endl;
      cout << "B= " << BA << endl;
      cout << "Chi2 = " << chi << endl; 
      cout << "----------------------------------------------------------" <<endl;
      cout << endl;
      cout << endl;

      results << "-------------------------NEW BIN-----------------------" << endl;
      //results << "bin definition:  "<< cutstring << endl;
      results << "Pt in " << jpsi_pt_bins[j] << ","<< jpsi_pt_bins[j+1] << endl;
      results << "Y in " << jpsi_y_bins[k] << "," << jpsi_y_bins[k+1]  << endl;
      results << "Pt bin center " << JPT << "   Y bin center " << JY << endl;
      results << "Name in plot  " << hst << endl;
      results << "-------------------------------------------------------" <<endl;
      results << "Fit : " << NSigA << " +/- " << errSigA << endl;
      results << "Resolution : " << resolA*1000. << " MeV  +/- "<<errresolA*1000  << endl; 
      results << "Mean : " << meanA << " GeV  +/- "<<errmeanA  << endl;
      results << "S/B= " << SNA << endl; 
      results << "Chi2 = " << chi << endl;
      results << "-------------------------------------------------------" <<endl;
      results << "----------------------------------J/psi-------------------------------" << endl;
      rfr->printMultiline(results,npar);
      results << "-------------------------------------------------------" <<endl;
      results << endl;
      results << endl;         
      delete ws;
      }
    }
    results.close();  
    out.cd();
    Hresults[k]->Write();
    HRMS[k]->Write();

    TCanvas cres;
    cres.cd();
    cres.SetLogy(1);
    Hresults[k]->SetMarkerSize(1.2);
    Hresults[k]->Draw("TEXT");
    sprintf(name,"InBins/ResultsJpsi_Y%i.png",k+1);
    cres.SaveAs(name);
    HRMS[k]->SetMarkerSize(1.2);
    HRMS[k]->Draw("TEXT");
    sprintf(name,"InBins/PtMeanRMS_Y%i.png",k+1);
    cres.SaveAs(name);
    delete Hresults[k];
    delete HRMS[k];
    out.Close();
  }

  //out.Close();
  return 1;
}
