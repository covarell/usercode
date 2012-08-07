// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2F.h>

void fitParsPt(int LHCsqrts = 7, bool isVBFsignal = false) {

  static const unsigned int Nmasspoints = 5;

  string masspointsS[Nmasspoints] = {"125","150","200","300","400"};
  if (isVBFsignal && LHCsqrts == 7) masspointsS[0] = "120";  

  double masspoints[Nmasspoints] = {0.,0.,0.,0.,0.};
  double masspointserr[Nmasspoints] = {0.,0.,0.,0.,0.};

  double bb[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Ti[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double enne[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double bbs[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tis[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double ennes[Nmasspoints] = {0.,0.,0.,0.,0.}; 

  double bberr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tierr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double enneerr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double bbserr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double Tiserr[Nmasspoints] = {0.,0.,0.,0.,0.}; 
  double enneserr[Nmasspoints] = {0.,0.,0.,0.,0.}; 

  ofstream *allParamsSig;
  ofstream *allParamsBkg;
  char fileName[200];

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    masspoints[i] = (double)atof(masspointsS[i].c_str());
  }
  
  if (isVBFsignal) sprintf(fileName,"allParamsVbf_%dTeV.txt",LHCsqrts);
  else sprintf(fileName,"allParamsSig_%dTeV.txt",LHCsqrts);
  allParamsSig = new ofstream(fileName);
  sprintf(fileName,"allParamsBkg_%dTeV.txt",LHCsqrts);
  allParamsBkg = new ofstream(fileName);

  char thePar[10];
  char equalS[1];
  char dashS[1];
  char pmS[3];
  char theLimit1[10];
  char theLimit2[10];
  float fitted, error;

  // Inizia
  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    if (isVBFsignal) sprintf(fileName,"text/paramsVbf_%dGeV_%dTeV_all.txt",int(masspoints[i]),LHCsqrts);
    else sprintf(fileName,"text/paramsSig_%dGeV_%dTeV_all.txt",int(masspoints[i]),LHCsqrts);
    ifstream theFile(fileName);
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      // limit the error to meaningful values
      if (fabs(error/fitted) > 0.3) error = 0.3*fitted;
      if (fabs(error/fitted) < 0.01) error = 0.01*fitted;
      if (!strcmp(thePar,"bbs")) {bbs[i] = fitted; bbserr[i] = error;}
      if (!strcmp(thePar,"Ts")) {Tis[i] = fitted; Tiserr[i] = error;}
      if (!strcmp(thePar,"ns")) {ennes[i] = fitted; enneserr[i] = error;}
    }

    theFile.close();

  }

  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    sprintf(fileName,"text/paramsBkg_%dGeV_%dTeV_all.txt",int(masspoints[i]),LHCsqrts);
    ifstream theFile(fileName);
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >>  dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      // limit the error to meaningful values
      if (fabs(error/fitted) > 0.3) error = 0.3*fitted;
      if (fabs(error/fitted) < 0.01) error = 0.01*fitted;
      if (!strcmp(thePar,"bb")) {bb[i] = fitted; bberr[i] = error;}
      if (!strcmp(thePar,"T")) {Ti[i] = fitted; Tierr[i] = error;}
      if (!strcmp(thePar,"n")) {enne[i] = fitted; enneerr[i] = error;}
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

  TGraphErrors *bbsfit = new TGraphErrors(Nmasspoints,masspoints,bbs,masspointserr,bbserr);
  
  bbsfit->SetTitle("bbs");
  bbsfit->SetMarkerStyle(20);
  bbsfit->SetMarkerColor(kRed);
  bbsfit->Draw("AP"); 
  bbsfit->Fit("mypol2","","PSAME"); 
  for (unsigned int j = 0; j < 3; j++) {
    *allParamsSig << "bbs" << j << " = " << mypol2->GetParameter(j) << " C \n";
  }

  Signal.cd(2);
  
  if (!isVBFsignal) {
    TGraphErrors *Tisfit = new TGraphErrors(Nmasspoints,masspoints,Tis,masspointserr,Tiserr);
    
    Tisfit->SetTitle("Tis");
    Tisfit->SetMarkerStyle(20);
    Tisfit->SetMarkerColor(kRed);
    Tisfit->Draw("AP"); 
    Tisfit->Fit("mypol2","","PSAME");
    for (unsigned int j = 0; j < 3; j++) {
      *allParamsSig << "Ts" << j << " = " << mypol2->GetParameter(j) << " C \n";
    }
    
    *allParamsSig << "ns0 = 0.733 C \n";
    *allParamsSig << "ns1 = 0. C \n";
    *allParamsSig << "ns2 = 0. C \n";
    *allParamsSig << "ms = 1066.2 C \n";
    *allParamsSig << "ndues = 0.95 C \n";

  } else {

    TGraphErrors *ennesfit = new TGraphErrors(Nmasspoints,masspoints,ennes,masspointserr,enneserr);
    
    ennesfit->SetTitle("ennes");
    ennesfit->SetMarkerStyle(20);
    ennesfit->SetMarkerColor(kRed);
    ennesfit->Draw("AP"); 
    ennesfit->Fit("mypol2","","PSAME");
    for (unsigned int j = 0; j < 3; j++) {
      *allParamsSig << "ns" << j << " = " << mypol2->GetParameter(j) << " C \n";
    }
    
    *allParamsSig << "Ts0 = 0.0064 C \n";
    *allParamsSig << "Ts1 = 0. C \n";
    *allParamsSig << "Ts2 = 0. C \n";
    *allParamsSig << "ms = 263.9 C \n";
    *allParamsSig << "ndues = 3.948 C \n";
  }

  Signal.SaveAs("fitPar_Signal.ps");

  allParamsSig->close();
  allParamsBkg->close();

  return;
}
