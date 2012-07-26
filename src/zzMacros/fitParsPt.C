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
#include <TH2F.h>

void fitParsPt() {

  static const unsigned int Nmasspoints = 5;

  string masspointsS[Nmasspoints] = {"125","150","200","300","400"};

  double masspoints[Nmasspoints] = {0.,0.,0.,0.,0.};
  double masspointserr[Nmasspoints] = {0.,0.,0.,0.,0.};

  double bb[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Ti[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double enne[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double bbs[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tis[Nmasspoints] = {0.,0.,0.,0.,0.}; 

  double bberr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tierr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double enneerr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double bbserr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tiserr[Nmasspoints] = {0.,0.,0.,0.,0.}; 

  ofstream *allParamsSig;
  ofstream *allParamsBkg;
  char fileout[200];

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    masspoints[i] = (double)atof(masspointsS[i].c_str());
  }
  allParamsSig = new ofstream("allParamsSig.txt");
  allParamsBkg = new ofstream("allParamsBkg.txt");

  string fileName;  
  char thePar[10];
  char equalS[1];
  char dashS[1];
  char pmS[3];
  char theLimit1[10];
  char theLimit2[10];
  float fitted, error;

  // Inizia
  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    fileName = "text/paramsSig_" + masspointsS[i] + "GeV_all.txt";
    ifstream theFile(fileName.c_str());
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      if (!strcmp(thePar,"bbs")) {bbs[i] = fitted; bbserr[i] = error;}
      if (!strcmp(thePar,"Ts")) {Tis[i] = fitted; Tiserr[i] = error;}
    }

    theFile.close();

  }

  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    fileName = "text/paramsBkg_" + masspointsS[i] + "GeV_all.txt";
    ifstream theFile(fileName.c_str());
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >>  dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      if (!strcmp(thePar,"bb")) {bb[i] = fitted; bberr[i] = error/10.;}
      if (!strcmp(thePar,"T")) {Ti[i] = fitted; Tierr[i] = error/10.;}
      if (!strcmp(thePar,"n")) {enne[i] = fitted; enneerr[i] = error/10.;}
    }

    theFile.close();

  }
  
  /* TF1 *myfunc1 = new TF1("myfunc1",slope1,200,1100,3);
  myfunc1->SetParameter(0,800);
  TF1 *myfunc2 = new TF1("myfunc2",slope121,200,1100,7);
  myfunc2->FixParameter(0,500);
  myfunc2->FixParameter(1,800);*/
  TF1 *mypol2 = new TF1("mypol2","pol2");
  TF1 *mypol3 = new TF1("mypol3","pol3");
  // TF1 *mypol4 = new TF1("mypol4","pol4");

 
  // A - Background
  TCanvas Background("Background","Background",10,10,1000,600);
  Background.Divide(2,2);
  Background.cd(1);
  
  TGraphErrors *bbfit = new TGraphErrors(Nmasspoints,masspoints,bb,masspointserr,bberr);

  bbfit->SetTitle("bb");
  bbfit->SetMarkerStyle(20);
  bbfit->SetMarkerColor(kRed);
  bbfit->Draw("AP"); 
  bbfit->Fit("mypol2","","PSAME"); 
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsBkg << "bb" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }

  Background.cd(2);
  
  TGraphErrors *Tifit = new TGraphErrors(Nmasspoints,masspoints,Ti,masspointserr,Tierr);

  Tifit->SetTitle("Ti");
  Tifit->SetMarkerStyle(20);
  Tifit->SetMarkerColor(kRed);
  Tifit->Draw("AP"); 
  Tifit->Fit("mypol2","","PSAME");
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsBkg << "T" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }

  Background.cd(3);
  
  TGraphErrors *ennefit = new TGraphErrors(Nmasspoints,masspoints,enne,masspointserr,enneerr);

  ennefit->SetTitle("enne");
  ennefit->SetMarkerStyle(20);
  ennefit->SetMarkerColor(kRed);
  ennefit->Draw("AP");
  ennefit->Fit("mypol2","","PSAME");
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsBkg << "n" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }
  *allParamsBkg << "m = 54.47 C \n";
  *allParamsBkg << "ndue = 1.73 C \n";

  Background.SaveAs("fitPar_Background.ps");

  // B - Signal
  TCanvas Signal("Signal","Signal",10,10,1000,300);
  Signal.Divide(2,1);

  Signal.cd(1);
  
  TGraphErrors *Tisfit = new TGraphErrors(Nmasspoints,masspoints,Tis,masspointserr,Tiserr);

  Tisfit->SetTitle("Tis");
  Tisfit->SetMarkerStyle(20);
  Tisfit->SetMarkerColor(kRed);
  Tisfit->Draw("AP"); 
  Tisfit->Fit("mypol2","","PSAME");
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsSig << "Ts" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }

  Signal.cd(2);
  TGraphErrors *bbsfit = new TGraphErrors(Nmasspoints,masspoints,bbs,masspointserr,bbserr);

  bbsfit->SetTitle("bbs");
  bbsfit->SetMarkerStyle(20);
  bbsfit->SetMarkerColor(kRed);
  bbsfit->Draw("AP"); 
  bbsfit->Fit("mypol2","","PSAME"); 
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsSig << "bbs" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }
  *allParamsSig << "ns0 = 0.733 C \n";
  *allParamsSig << "ns1 = 0. C \n";
  *allParamsSig << "ns2 = 0. C \n";
  *allParamsSig << "ms = 1066.2 C \n";
  *allParamsSig << "ndues = 0.95 C \n";

  Signal.SaveAs("fitPar_Signal.ps");

  allParamsSig->close();
  allParamsBkg->close();

  return;
}
