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

Double_t slope1(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= par[0]) fitval = par[1];
  else fitval = par[2]*x[0] + par[1] - par[2]*par[0];
  return fitval;
}

Double_t slope121(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= par[0]) fitval = par[5]*x[0] - par[5]*par[0] + par[2]*par[0]*par[0] + par[3]*par[0] - par[4];
  else if (x[0] >= par[1]) fitval = par[6]*x[0] - par[6]*par[1] + par[2]*par[1]*par[1] + par[3]*par[1] - par[4];
  else fitval = par[2]*x[0]*x[0] + par[3]*x[0] - par[4];
  return fitval;
}

void fitAccPars() {

  static const unsigned int Nmasspoints = 6;

  string masspointsS[Nmasspoints] = {"300","500","600","700","800","1000"};

  double masspoints[Nmasspoints] = {0.,0.,0.,0.,0.,0.};
  double masspointserr[Nmasspoints] = {0.,0.,0.,0.,0.,0.};

  double para2[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para4[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para6[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para8[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double acca2[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double a2[Nmasspoints] = {0.,0.,0.,0.,0.,0.};
  double cutOff[Nmasspoints] = {0.,0.,0.,0.,0.,0.};
  double g[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double a4[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double b2[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double b4[Nmasspoints] = {0.,0.,0.,0.,0.,0.};

  double para2err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para4err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para6err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double para8err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double acca2err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double a2err[Nmasspoints] = {0.,0.,0.,0.,0.,0.};
  double cutOfferr[Nmasspoints] = {0.,0.,0.,0.,0.,0.};
  double gerr[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double a4err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double b2err[Nmasspoints] = {0.,0.,0.,0.,0.,0.}; 
  double b4err[Nmasspoints] = {0.,0.,0.,0.,0.,0.};

  ofstream *allParams;
  ofstream *of[Nmasspoints];
  char fileout[200];

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    masspoints[i] = (double)atof(masspointsS[i].c_str());
    sprintf(fileout,"acc_parsEstim%s.txt",masspointsS[i].c_str());
    of[i] = new ofstream(fileout);
  }
  allParams = new ofstream("allParams.txt");

  string fileName;  
  char thePar[10];
  float fitted, error;

  // Inizia
  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    fileName = "acc_parsFinal" + masspointsS[i] + ".txt";
    ifstream theFile(fileName.c_str());
	
    while (theFile >> thePar >> fitted >> error) {
      cout << thePar << " " << fitted << " " << error << endl;
      if (error < 0.001) error = 0.01; 
      if (!strcmp(thePar,"para2")) {para2[i] = fitted; para2err[i] = error;}
      if (!strcmp(thePar,"para4")) {para4[i] = fitted; para4err[i] = error;}
      if (!strcmp(thePar,"para6")) {para6[i] = fitted; para6err[i] = error;}
      if (!strcmp(thePar,"para8")) {para8[i] = fitted; para8err[i] = error;}
      if (!strcmp(thePar,"a2")) {a2[i] = fitted; a2err[i] = error;}
      if (!strcmp(thePar,"a4")) {a4[i] = fitted; a4err[i] = error;}
      if (!strcmp(thePar,"g")) {g[i] = fitted; gerr[i] = error;}
      if (!strcmp(thePar,"cutOff")) {cutOff[i] = fitted; cutOfferr[i] = error;}
      if (!strcmp(thePar,"b2")) {b2[i] = fitted; b2err[i] = error;}
      if (!strcmp(thePar,"b4")) {b4[i] = fitted; b4err[i] = error;}
      if (!strcmp(thePar,"acca2")) {acca2[i] = fitted; acca2err[i] = error;}
    }

    theFile.close();

  }
  
  TF1 *myfunc1 = new TF1("myfunc1",slope1,200,1100,3);
  myfunc1->SetParameter(0,800);
  TF1 *myfunc2 = new TF1("myfunc2",slope121,200,1100,7);
  myfunc2->FixParameter(0,500);
  myfunc2->FixParameter(1,800);
  // TF1 *mypol2 = new TF1("mypol2","pol2");
  TF1 *mypol3 = new TF1("mypol3","pol3");
  // TF1 *mypol4 = new TF1("mypol4","pol4");

 
  // A - Cos(theta*)
  TCanvas cosThetaStar("cosThetaStar","cosThetaStar",10,10,1000,600);
  cosThetaStar.Divide(2,2);
  cosThetaStar.cd(1);
  
  TGraphErrors *para2fit = new TGraphErrors(Nmasspoints,masspoints,para2,masspointserr,para2err);

  para2fit->SetTitle("para2");
  para2fit->SetMarkerStyle(20);
  para2fit->SetMarkerColor(kRed);
  para2fit->Draw("AP"); 
  para2fit->Fit("myfunc2","","PSAME"); 
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "[Acceptance] \n";
    *of[i] << "para2 = " << myfunc2->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 5; j++) {
    *allParams << "para2 " << j << " " << myfunc2->GetParameter(j+2) << "\n";
  }

  cosThetaStar.cd(2);
  
  TGraphErrors *para4fit = new TGraphErrors(Nmasspoints,masspoints,para4,masspointserr,para4err);

  para4fit->SetTitle("para4");
  para4fit->SetMarkerStyle(20);
  para4fit->SetMarkerColor(kRed);
  para4fit->Draw("AP"); 
  para4fit->Fit("myfunc2","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "para4 = " << myfunc2->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 5; j++) {
    *allParams << "para4 " << j << " " << myfunc2->GetParameter(j+2) << "\n";
  }

  cosThetaStar.cd(3);
  
  TGraphErrors *para6fit = new TGraphErrors(Nmasspoints,masspoints,para6,masspointserr,para6err);

  para6fit->SetTitle("para6");
  para6fit->SetMarkerStyle(20);
  para6fit->SetMarkerColor(kRed);
  para6fit->Draw("AP");
  para6fit->Fit("myfunc1","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "para6 = " << myfunc1->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 2; j++) {
    *allParams << "para6 " << j << " " << myfunc1->GetParameter(j+1) << "\n";
  }
 
  cosThetaStar.cd(4);
  
  /* TGraphErrors *para8fit = new TGraphErrors(Nmasspoints,masspoints,para8,masspointserr,para8err);

  para8fit->SetTitle("para8");
  para8fit->SetMarkerStyle(20);
  para8fit->SetMarkerColor(kRed);
  para8fit->Draw("AP");
  para8fit->Fit("myfunc2","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "para8 = " << myfunc2->Eval(masspoints[i]) << " C \n";
  }  */

  cosThetaStar.SaveAs("fitPar_cosThetaStar.ps");

  // B - Cos(theta1)
  TCanvas cosTheta1("cosTheta1","cosTheta1",10,10,500,600);
  cosTheta1.Divide(1,2);
  cosTheta1.cd(1);
  
  TGraphErrors *b2fit = new TGraphErrors(Nmasspoints,masspoints,b2,masspointserr,b2err);

  b2fit->SetTitle("b2");
  b2fit->SetMarkerStyle(20);
  b2fit->SetMarkerColor(kRed);
  b2fit->Draw("AP"); 
  b2fit->Fit("mypol3","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "b2 = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "b2 " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta1.cd(2);
  
  TGraphErrors *b4fit = new TGraphErrors(Nmasspoints,masspoints,b4,masspointserr,b4err);

  b4fit->SetTitle("b4");
  b4fit->SetMarkerStyle(20);
  b4fit->SetMarkerColor(kRed);
  b4fit->Draw("AP"); 
  b4fit->Fit("mypol3","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "b4 = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "b4 " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta1.SaveAs("fitPar_cosTheta1.ps");

  // C - Cos(theta2)
  TCanvas cosTheta2("cosTheta2","cosTheta2",10,10,1000,600);
  cosTheta2.Divide(2,2);
  cosTheta2.cd(1);
  
  TGraphErrors *a2fit = new TGraphErrors(Nmasspoints,masspoints,a2,masspointserr,a2err);

  a2fit->SetTitle("a2");
  a2fit->SetMarkerStyle(20);
  a2fit->SetMarkerColor(kRed);
  a2fit->Draw("AP"); 
  a2fit->Fit("mypol3","","PSAME"); 
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "a2 = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "a2 " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta2.cd(2);
  
  TGraphErrors *a4fit = new TGraphErrors(Nmasspoints,masspoints,a4,masspointserr,a4err);

  a4fit->SetTitle("a4");
  a4fit->SetMarkerStyle(20);
  a4fit->SetMarkerColor(kRed);
  a4fit->Draw("AP"); 
  a4fit->Fit("mypol3","","PSAME"); 
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "a4 = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "a4 " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta2.cd(3);
  
  TGraphErrors *gfit = new TGraphErrors(Nmasspoints,masspoints,g,masspointserr,gerr);

  gfit->SetTitle("g");
  gfit->SetMarkerStyle(20);
  gfit->SetMarkerColor(kRed);
  gfit->Draw("AP");
  gfit->Fit("mypol3","","PSAME"); 
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "g = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "g " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta2.cd(4);
  
  TGraphErrors *cutOfffit = new TGraphErrors(Nmasspoints,masspoints,cutOff,masspointserr,cutOfferr);

  cutOfffit->SetTitle("cutOff");
  cutOfffit->SetMarkerStyle(20);
  cutOfffit->SetMarkerColor(kRed);
  cutOfffit->Draw("AP");
  cutOfffit->Fit("mypol3","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "cutOff = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "cutOff " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  cosTheta2.SaveAs("fitPar_cosTheta2.ps");

  // D - phi1
  TCanvas phi1("phi1","phi1",10,10,500,300);
  phi1.cd();
  
  TGraphErrors *acca2fit = new TGraphErrors(Nmasspoints,masspoints,acca2,masspointserr,acca2err);

  acca2fit->SetTitle("acca2");
  acca2fit->SetMarkerStyle(20);
  acca2fit->SetMarkerColor(kRed);
  acca2fit->Draw("AP"); 
  acca2fit->Fit("mypol3","","PSAME");
  for (unsigned int i = 0; i < Nmasspoints; i++) {
    *of[i] << "acca2 = " << mypol3->Eval(masspoints[i]) << " C \n";
  }
  for (unsigned int j = 0; j < 4; j++) {
    *allParams << "acca2 " << j << " " << mypol3->GetParameter(j) << "\n";
  }

  phi1.SaveAs("fitPar_phi1.ps");

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    of[i]->close();
  }
  allParams->close();

  return;
}
