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

Double_t zzFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 182.) fitval = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  else fitval = par[3]*x[0] + par[4]*x[0]*x[0] - par[3]*182. - par[4]*33124. + par[0] + par[1]*182. + par[2]*33124.;
  return fitval;
}

void fitParsPtCJLST(int LHCsqrts = 7, int whichtype = 1) {

// whichtype
// 0 - gg Signal
// 1 - VBF Signal
// 2 - ZZ
// 3 - ZX
// 4 - ggZZ
// 5 - WH
// 6 - ZH
// 7 - ttH  

  cout << "FIRST RUN source changeCinZero.csh!" << endl;

  static const unsigned int Nmasspoints = 9;
  string masspointsS[Nmasspoints] = {"115","120","125","130","140","200","400","700","1000"};

  // static const unsigned int Nmasspoints = 7;
  // string masspointsS[Nmasspoints] = {"120","130","140","200","400","700","1000"};
  string nameSample[8] = {"gg","vbf","zz","zx","ggzz","wh","zh","tth"};

  double masspoints[Nmasspoints] = {Nmasspoints*0.};
  double masspointserr[Nmasspoints] = {Nmasspoints*0.};

  double bb[Nmasspoints] = {Nmasspoints*0.}; 
  double T[Nmasspoints] = {Nmasspoints*0.}; 
  double n[Nmasspoints] = {Nmasspoints*0.}; 
  double bbdue[Nmasspoints] = {Nmasspoints*0.}; 
  double ndue[Nmasspoints] = {Nmasspoints*0.}; 
  double m[Nmasspoints] = {Nmasspoints*0.}; 
  double fexp[Nmasspoints] = {Nmasspoints*0.}; 

  double bberr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double Terr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double nerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double bbdueerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double ndueerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double merr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double fexperr[Nmasspoints] = {Nmasspoints*0.000000001}; 

  ofstream *allParams;
  char fileName[200];

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    masspoints[i] = (double)atof(masspointsS[i].c_str());
  }
  
  sprintf(fileName,"text/allParams_%s_%dTeV.txt",nameSample[whichtype].c_str(),LHCsqrts);
  allParams = new ofstream(fileName);

  char thePar[10];
  char equalS[1];
  char dashS[1];
  char pmS[3];
  char theLimit1[10];
  char theLimit2[10];
  float fitted, error;

  // Inizia
  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_Default.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts);
    ifstream theFile(fileName);
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      // limit the error to meaningful values
      // if (fabs(error/fitted) > 0.3) error = 0.3*fitted;
      // if (fabs(error/fitted) < 0.01) error = 0.01*fitted;
      if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) {bb[i] = fitted; bberr[i] = error;}
      if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) {T[i] = fitted; Terr[i] = error;}
      if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) {bbdue[i] = fitted; bbdueerr[i] = error;}
      if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) {fexp[i] = fitted; fexperr[i] = error;}
      if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) {n[i] = fitted; nerr[i] = error;}
      if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) {m[i] = fitted; merr[i] = error;}
      if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) {ndue[i] = fitted; ndueerr[i] = error;}
    }

    theFile.close();

  }

  TF1 *mypol0 = new TF1("mypol0","pol0");
  TF1 *mypol2 = new TF1("mypol2","pol2");
  TF1 *mypol3 = new TF1("mypol3","pol3");
  TF1 *myZZfunc = new TF1("myZZfunc",zzFunc,100.,1600.,5);
  // TF1 *mypol4 = new TF1("mypol4","pol4");
 
  TCanvas Background("Background","Background",10,10,600,1000);
  Background.Divide(2,4);

  Background.cd(1);
  
  int nPointsRemoved = 0;

  TGraphErrors *bbfit = new TGraphErrors(Nmasspoints,masspoints,bb,masspointserr,bberr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(bb[i]) < 0.0000001) {
      cout << "Point of mass " << masspointsS[i] << " is not there, removed" << endl; 
      bbfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  bbfit->SetTitle("bb");
  bbfit->SetMarkerStyle(20);
  bbfit->SetMarkerColor(kRed);
  bbfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (bberr[1] > 0.0000001 && bberr[Nmasspoints-1] > 0.0000001) {  
    bbfit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "bb [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    bbfit->Fit("mypol0","","PSAME"); 
    *allParams << "bb [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(2);
  nPointsRemoved = 0;

  TGraphErrors *Tfit = new TGraphErrors(Nmasspoints,masspoints,T,masspointserr,Terr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(T[i]) < 0.0000001) {
      Tfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;     
    }
  }
  
  Tfit->SetTitle("T");
  Tfit->SetMarkerStyle(20);
  Tfit->SetMarkerColor(kRed);
  Tfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (Terr[1] > 0.0000001 && Terr[Nmasspoints-1] > 0.0000001) { 
    Tfit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "T [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    Tfit->Fit("mypol0","","PSAME"); 
    *allParams << "T [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(3);
  nPointsRemoved = 0;

  TGraphErrors *bbduefit = new TGraphErrors(Nmasspoints,masspoints,bbdue,masspointserr,bbdueerr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(bbdue[i]) < 0.0000001) {
      bbduefit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  bbduefit->SetTitle("bb2");
  bbduefit->SetMarkerStyle(20);
  bbduefit->SetMarkerColor(kRed);
  bbduefit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (bbdueerr[1] > 0.0000001 && bbdueerr[Nmasspoints-1] > 0.0000001) { 
    bbduefit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "bbdue [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    bbduefit->Fit("mypol0","","PSAME"); 
    *allParams << "bbdue [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(4);
  nPointsRemoved = 0;

  TGraphErrors *fexpfit = new TGraphErrors(Nmasspoints,masspoints,fexp,masspointserr,fexperr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(fexp[i]) < 0.0000001) {
      fexpfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  fexpfit->SetTitle("fexp");
  fexpfit->SetMarkerStyle(20);
  fexpfit->SetMarkerColor(kRed);
  fexpfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (fexperr[1] > 0.0000001 && fexperr[Nmasspoints-1] > 0.0000001) { 
    fexpfit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "fexp [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    fexpfit->Fit("mypol0","","PSAME"); 
    *allParams << "fexp [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(5);
  nPointsRemoved = 0;

  TGraphErrors *nfit = new TGraphErrors(Nmasspoints,masspoints,n,masspointserr,nerr);
  nfit->Print();
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(n[i]) < 0.0000001) {
      nfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  nfit->SetTitle("n");
  nfit->SetMarkerStyle(20);
  nfit->SetMarkerColor(kRed);
  nfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (nerr[1] > 0.0000001 && nerr[Nmasspoints-1] > 0.0000001) { 
    nfit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "n [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    nfit->Fit("mypol0","","PSAME"); 
    *allParams << "n [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(6);
  nPointsRemoved = 0;

  TGraphErrors *mfit = new TGraphErrors(Nmasspoints,masspoints,m,masspointserr,merr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(m[i]) < 0.0000001) {
      mfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  mfit->SetTitle("m");
  mfit->SetMarkerStyle(20);
  mfit->SetMarkerColor(kRed);
  mfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (merr[1] > 0.0000001 && merr[Nmasspoints-1] > 0.0000001) { 
    mfit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "m [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    mfit->Fit("mypol0","","PSAME"); 
    *allParams << "m [0] = " << mypol0->GetParameter(0) << "\n";
  }

  Background.cd(7);
  nPointsRemoved = 0;

  TGraphErrors *nduefit = new TGraphErrors(Nmasspoints,masspoints,ndue,masspointserr,ndueerr);
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(ndue[i]) < 0.0000001) {
      nduefit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  nduefit->SetTitle("n2");
  nduefit->SetMarkerStyle(20);
  nduefit->SetMarkerColor(kRed);
  nduefit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (ndueerr[1] > 0.0000001 && ndueerr[Nmasspoints-1] > 0.0000001) { 
    nduefit->Fit("myZZfunc","","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << "ndue [" << j << "] = " << myZZfunc->GetParameter(j) << "\n";
    }
  } else {
    nduefit->Fit("mypol0","","PSAME"); 
    *allParams << "ndue [0] = " << mypol0->GetParameter(0) << "\n";
  }

  sprintf(fileName,"figs/fitParCJLST_%s_%dTeV.pdf",nameSample[whichtype].c_str(),LHCsqrts);
  Background.SaveAs(fileName);
  allParams->close();

  return;
}
