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
  ws->factory("Polynomial::CPolFunct(JpsiMass,{CoefPol1[-1.0,-300.,1000.]})");
  ws->factory("Polynomial::CPolFunct2(JpsiMass,{CoefPol1,CoefPol2[0.5,0.,1000]})");
  ws->factory("Exponential::expFunct(JpsiMass,coefExp[-50,-80.,10])");

  return;
}

void defineSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:
  cout << "defining signal" << endl;
  //Normal Gaussians
  ws->factory("Gaussian::signalG1(JpsiMass,meanSig1[3.1,3.05,3.15],sigmaSig1[0.03,0.008,0.08])");
  ws->factory("Gaussian::signalG2(JpsiMass,meanSig1,sigmaSig2[0.02,0.008,0.1])");
  // ws->factory("Gaussian::signalG2(JpsiMass,meanSig1,sigmaSig1");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(JpsiMass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(JpsiMass,meanSig1,sigmaSig1,alpha[0.5,.1,3.],enne[8.,2.,10.])");

 ws->factory("CBShape::sigCB2(JpsiMass,meanSig1,sigmaSigCB2[0.04,0.,4.],alpha1[0.5,0.,4.],enne1[5.,1.,50.])");
  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  ws->factory("SUM::sigPDF(coeffGauss[0.01,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystallBall
  ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");


  //Sum of Two CBs
  ws->factory("SUM::sigDoubleCB(coeffGauss*sigCB,sigCB2)");

  return;
}

void prefitSideband(RooWorkspace *ws, const int DataCat)
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

  // ws->var("coefExp")->setConstant(kTRUE);
  // ws->var("alpha")->setConstant(kTRUE);
  // ws->var("enne")->setConstant(kTRUE);
  // ws->var("CoefPol1")->setConstant(kTRUE);
  // ws->var("CoefPol2")->setConstant(kTRUE);
  // ws->var("sigmaSig2")->setConstant(kTRUE);
  // ws->var("meanSig1")->setConstant(kTRUE);
  
  return;
}

void setRanges(RooWorkspace *ws)
{
  cout<< "setting ranges" << endl;
  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.5;

  ws->var("JpsiMass")->setRange("all",JpsiMassMin,JpsiMassMax);
  ws->var("JpsiMass")->setRange("left1",2.75,2.9);
  ws->var("JpsiMass")->setRange("left",JpsiMassMin,2.95);
  ws->var("JpsiMass")->setRange("right",3.25,JpsiMassMax);
  ws->var("JpsiMass")->setRange("right1",3.4,JpsiMassMax);
  ws->var("JpsiMass")->setRange("peak",2.95,3.15);

  ws->cat("JpsiType")->setRange("glbglb","GG");
  ws->cat("JpsiType")->setRange("glbtrk","GT");
  ws->cat("JpsiType")->setRange("trktrk","TT");

  return;
}

double drawResults(RooWorkspace *ws,const int npar, char *hst, TH2F * hrms, int bin)

{
  cout << "drawing results" << endl;
 
  RooDataSet *data = (RooDataSet*)ws->data("data");
  RooAbsPdf *totPDF = ws->pdf("totPDF");

  char reducestr[200];
  const double sigma = ws->var("sigmaSig1")->getVal();
  const double mean = ws->var("meanSig1")->getVal();
  double up=mean+(2.5*sigma);
  double down=mean-(2.5*sigma);

  RooPlot *mframe = ws->var("JpsiMass")->frame(2.6,3.5);
  RooPlot *ptframe = ws->var("JpsiPt")->frame();
  RooPlot *etaframe = ws->var("JpsiEta")->frame();


   mframe->SetTitle("");
   
   //Plotting Data    
   data->plotOn(mframe,DataError(RooAbsData::SumW2),Binning(90),MarkerStyle(7),MarkerSize(.5));
   data->statOn(ptframe);
   data->plotOn(ptframe,DataError(RooAbsData::SumW2),Cut("JpsiMass>= 3.0 && JpsiMass <= 3.2"),Binning(100));
   data->plotOn(etaframe,DataError(RooAbsData::SumW2),Cut("JpsiMass>= 3.0 && JpsiMass <= 3.2"),Binning(100));
   RooDataSet* datar=(RooDataSet*)data->reduce("JpsiMass>= 3.0 && JpsiMass <= 3.2");
   RooRealVar *meanv= datar->meanVar(*ws->var("JpsiPt"));
   RooRealVar *rmsv= datar->rmsVar(*ws->var("JpsiPt"));
   hrms->SetBinContent(bin,meanv->getVal());
   hrms->SetBinError(bin,rmsv->getVal());
  
   //Plotting PDF
   totPDF->plotOn(mframe,LineColor(kBlack),LineWidth(1.0),Normalization(1.0,RooAbsReal::RelativeExpected));
   double chi=mframe->chiSquare(npar);

   totPDF->plotOn(mframe,Components("expFunct"),LineColor(kBlue),LineWidth(1.0),LineStyle(kDashed),Normalization(1.0,RooAbsReal::RelativeExpected));

  
   mframe->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
 

   TCanvas c1;

   c1.cd();
   
   TLatex *t = new TLatex();
   t->SetNDC();
   t->SetTextAlign(13);
   t->SetTextFont(63);
   t->SetTextSizePixels(25);
   
   char chiSq[40];
   sprintf(chiSq,"#chi^{2}= %5.2f",chi);
   mframe->Draw();
   t->DrawLatex(0.2,0.85,chiSq);

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

void printResults(RooWorkspace *ws, double &Nsig, double &errSig, double &resol,double &Mean, double &errMean, double &ErrResol, double &SN, double &Bkg, TH2F * h,int bin)
{
  cout << "printing " << endl;
  Nsig   = ws->var("NSig")->getVal();
  double Nbkg = ws->var("NBkg")->getVal();
  errSig = ws->var("NSig")->getError();
  double coeffGauss = ws->var("coeffGauss")->getVal();
  double sigmaSig1 = ws->var("sigmaSig1")->getVal();
  ErrResol=ws->var("sigmaSig1")->getError();
  Mean = ws->var("meanSig1")->getVal();
  errMean=ws->var("meanSig1")->getError();
  double sigmaSig2 = ws->var("sigmaSig2")->getVal();
  resol = sigmaSig1;

 ws->var("JpsiMass")->setRange("integral",Mean-2.5*resol,Mean+2.5*resol);

 RooAbsReal* N= ws->pdf("expFunct")->createIntegral(RooArgSet(*ws->var("JpsiMass")),NormSet(RooArgSet(*ws->var("JpsiMass"))),Range("integral"));
 RooAbsReal* S= ws->pdf("sigCBGauss")->createIntegral(RooArgSet(*ws->var("JpsiMass")),NormSet(RooArgSet(*ws->var("JpsiMass"))),Range("integral"));
  
  Bkg=Nbkg*N->getVal();
						   
  SN=(Nsig*S->getVal())/(Nbkg*N->getVal());
  h->SetBinContent(bin,Nsig);
  h->SetBinError(bin,errSig);
 
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
 
  Double_t jpsi_y_bins[] = {0., 1.2, 1.6, 2.4}; 

  Double_t jpsi_pt_bins0[] = {4.5, 5.5, 6.5, 8, 10, 12, 30};
  Double_t jpsi_pt_bins1[] = {2, 3.5, 4.5, 5.5, 6.5, 8, 10, 30};
  Double_t jpsi_pt_bins2[] = {0.,0.5,0.75,1,1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3, 3.25, 3.5, 4, 4.5, 5.5, 6.5, 8, 10,30};

  int sY=sizeof(jpsi_y_bins)/sizeof(double);
  int sp0=sizeof(jpsi_pt_bins0)/sizeof(double);
  int sp1=sizeof(jpsi_pt_bins1)/sizeof(double);
  int sp2=sizeof(jpsi_pt_bins2)/sizeof(double);
 
  
  /* TH2F * Hresults= new TH2F("Yields","Yields",JPSI_Y_BINS,jpsi_y_bins,JPSI_PT_BINS,jpsi_pt_bins);
  TH2F * HRMS= new TH2F("PtMedioRms","PtMedioRms",JPSI_Y_BINS,jpsi_y_bins,JPSI_PT_BINS,jpsi_pt_bins);
  double JY,JPT;
  */
  char name[300];
  sprintf(name,"InBins/FitResults.txt"); 
  std::ofstream results;
  results.open(name);
 
  TFile fIn(filename);
  fIn.cd();
  RooDataSet *reddata2 = (RooDataSet*)fIn.Get("data");

  sprintf(name,"InBins/HistoResults.root");
  TFile out(name,"RECREATE");

  for (int k=0; k<sY-1; k++){

    int sizeP;
    Double_t * jpsi_pt_bins;
    if (k==0) {jpsi_pt_bins = jpsi_pt_bins0;sizeP=sp0;}
    if (k==1) {jpsi_pt_bins = jpsi_pt_bins1;sizeP=sp1;}
    if (k==2) {jpsi_pt_bins = jpsi_pt_bins2;sizeP=sp2;}

    sprintf(name,"YieldsY%i",k);   
    TH2F * Hresults= new TH2F(name,name,sY-1,jpsi_y_bins,sizeP-1,jpsi_pt_bins);
    sprintf(name,"PtMedioRMS%i",k);
    TH2F * HRMS= new TH2F(name,name,sY-1,jpsi_y_bins,sizeP-1,jpsi_pt_bins);
    double JY,JPT;
   

    for (int j=0; j<sizeP-1; j++){

      JY=(jpsi_y_bins[k]+jpsi_y_bins[k+1])/2;
      JPT=(jpsi_pt_bins[j]+jpsi_pt_bins[j+1])/2;

      char cutstring[300];
      char hst[30];
  
      sprintf(cutstring,"((JpsiEta > -%5.2f && JpsiEta <= -%5.2f) || (JpsiEta >= %5.2f && JpsiEta < %5.2f)) && (JpsiPt >= %5.2f && JpsiPt < %5.2f)",jpsi_y_bins[k+1],jpsi_y_bins[k],jpsi_y_bins[k],jpsi_y_bins[k+1],jpsi_pt_bins[j],jpsi_pt_bins[j+1]); 
      
      sprintf(hst,"PT%i_Y%i",j+1,k+1);
      //cout << " Cut on data: " << cutstring << endl;

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
  
      ws->factory("SUM::totPDF(NSig[5000.,4.,10000000.]*sigCBGauss,NBkg[2000.,4.,10000000.]*expFunct)");
 
      int bin=Hresults->FindBin(JY,JPT);
 
      if (sidebandPrefit) prefitSideband(ws,k);

      RooFitResult *rfr =ws->pdf("totPDF")->fitTo(*reddata,Extended(1),Save(1),Minos(0),NumCPU(6),SumW2Error(kTRUE));
     
      int npar=rfr->floatParsFinal().getSize();
      double NSigA, errSigA,resolA,meanA,errmeanA,errresolA,SNA,chi,BA;
      chi=drawResults(ws,npar,hst,HRMS,bin);
      printResults(ws,NSigA,errSigA,resolA,meanA,errmeanA,errresolA,SNA,BA,Hresults,bin);

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

    TCanvas cres;
    cres.cd();
    //cres.SetLogy(1);
    Hresults->Draw("TEXT");
    sprintf(name,"InBins/ResultsY%i.png",k);
    cres.SaveAs(name);
    HRMS->Draw("TEXT");
    sprintf(name,"InBins/PtMeanRMSY%i.png",k);
    cres.SaveAs(name);
    out.cd();
    Hresults->Write();
    HRMS->Write();
    delete HRMS;
    delete Hresults;
  }

  out.Close();
  return 1;
}
