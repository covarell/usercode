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

bool isEfficiency = true;

void do2DResultSystFineBins(string whichTypeJpsi = "BP", string fileNameBase = "results", string fileNameComp = "results2", int maxBetween = 1, string fileNameComp2 = "none", string fileNameComp3 = "none") {

  double ybinlimits[6] = {0.0,0.9,1.2,1.6,2.1,2.4};

  static const unsigned int nbinspt1 = 11;
  static const unsigned int nbinspt2 = 7;
  static const unsigned int nbinspt3 = 12;
  static const unsigned int nbinspt4 = 12;
  static const unsigned int nbinspt5 = 6;

  double pt1binlimits[nbinspt1+1] = {7.,8.,9.,10.,11.,12.,13.5,15.,18.,30.,45.,70.};
  double pt2binlimits[nbinspt2+1] = {6.5,8.,9.,10.,12.,15.,30.,45.};
  double pt3binlimits[nbinspt3+1] = {5.5,6.5,7.,7.5,8.,8.5,9.,10.,11.,12.,15.,30.,45.};
  double pt4binlimits[nbinspt4+1] = {5.5,6.5,7.,7.25,7.5,8.,8.5,9.,10.,11.,12.,15.,30.};
  double pt5binlimits[nbinspt5+1] = {5.5,8.,9.,10.,12.,15.,30.};
  
  TH2F histJpsi00_09("histJpsi00_09","histJpsi",5,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histJpsi09_12("histJpsi09_12","histJpsi",5,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histJpsi12_16("histJpsi12_16","histJpsi",5,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histJpsi16_21("histJpsi16_21","histJpsi",5,ybinlimits,nbinspt4,pt4binlimits);
  TH2F histJpsi21_24("histJpsi21_24","histJpsi",5,ybinlimits,nbinspt5,pt5binlimits);

  TH2F histRelVar00_09("histRelVar00_09","histJpsi - relative variation",5,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histRelVar09_12("histRelVar09_12","histJpsi - relative variation",5,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histRelVar12_16("histRelVar12_16","histJpsi - relative variation",5,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histRelVar16_21("histRelVar16_21","histJpsi - relative variation",5,ybinlimits,nbinspt4,pt4binlimits);
  TH2F histRelVar21_24("histRelVar21_24","histJpsi - relative variation",5,ybinlimits,nbinspt5,pt5binlimits);

  gROOT->ProcessLine(".! rm -f lista lista2 lista3");
  char theCommand[200];

  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  sprintf(theCommand,".! ls results/%s*.txt > lista2",fileNameComp.c_str());
  gROOT->ProcessLine(theCommand);  
  ifstream lista("lista");
  ifstream lista2("lista2");
  ifstream lista3, lista4;
  if (maxBetween >= 2) {
    sprintf(theCommand,".! ls results/%s*.txt > lista3",fileNameComp2.c_str());
    gROOT->ProcessLine(theCommand);  
    lista3.open("lista3");
  }
   if (maxBetween == 3) {
    sprintf(theCommand,".! ls results/%s*.txt > lista4",fileNameComp3.c_str());
    gROOT->ProcessLine(theCommand);  
    lista4.open("lista4");
  }

  char fileName[300];
  char fileName2[300];
  char fileName3[300];
  char fileName4[300];
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  char theType2[2];
  float trueMC, fitted, fitted2, fitted3, fitted4, error, error2; 
  int i = 0;       
  string cutstring;
  float chi2i = 0.;  
  float theMaximumPt = -1.;
  float theMinimumPt = 9000000.;
  int theBin = 0;
  // Inizia

  while (lista >> fileName) {

    // scroll lines in the same order
    lista2 >> fileName2;
    if (maxBetween >= 2) lista3 >> fileName3;
    if (maxBetween == 3) lista4 >> fileName4;

    cutstring = "results/" + fileNameBase + "_pT%f-%f_y%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
      ifstream theFile2(fileName2);
      ifstream theFile3, theFile4;
      if (maxBetween >= 2) theFile3.open(fileName3);
      if (maxBetween == 3) theFile4.open(fileName4);

      while (theFile >> theType >> trueMC >> fitted >> error) {
	
	// scroll lines in the same order
	theFile2 >> theType2 >> trueMC >> fitted2 >> error2;
	cout << "TEST : " << theType << " " << trueMC << " " << fitted << " " << fitted2 << " " << error << endl;
	if (maxBetween >= 2) theFile3 >> theType2 >> trueMC >> fitted3 >> error2;
	if (maxBetween == 3) theFile4 >> theType2 >> trueMC >> fitted4 >> error2;

	if (!strcmp(theType,whichTypeJpsi.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin + 0.01)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          float eff = fitted/fitted2;
	  if (eff < 0.0) eff = 0.0;   if (eff > 1.0) eff = 1.0;
	  float efferr = sqrt(eff*(1.0-eff)/fitted2);
	  float relvar = 100.0*fabs(fitted-fitted2)/fitted2;
	  if (maxBetween >= 2) {
	    float relvar2 = 100.0*fabs(fitted2-fitted3)/fitted2;
	    if (relvar2 > relvar) {
	      relvar = relvar2;   fitted = fitted3;
	    }
	  }
	  if (maxBetween == 3) {
	    float relvar3 = 100.0*fabs(fitted2-fitted4)/fitted2;
	    if (relvar3 > relvar) {
	      relvar = relvar3;   fitted = fitted4;
	    }
	  }
	  if (isEfficiency) {
	    fitted = 100.0*eff;    error = 100.0*efferr;
	  }
	  if (fabs(ybincenter-0.45) < 0.01 ) {
	    theBin = histJpsi00_09.FindBin(ybincenter,ptbincenter);
	    histJpsi00_09.SetBinContent(theBin,fitted);
	    histJpsi00_09.SetBinError(theBin,error);
	    histRelVar00_09.SetBinContent(theBin,relvar);
	  } else if (fabs(ybincenter-1.05) < 0.01 ) {
	    theBin = histJpsi09_12.FindBin(ybincenter,ptbincenter);
	    histJpsi09_12.SetBinContent(theBin,fitted);
	    histJpsi09_12.SetBinError(theBin,error); 
	    histRelVar09_12.SetBinContent(theBin,relvar);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histJpsi12_16.FindBin(ybincenter,ptbincenter);
	    histJpsi12_16.SetBinContent(theBin,fitted);
	    histJpsi12_16.SetBinError(theBin,error); 
	    histRelVar12_16.SetBinContent(theBin,relvar);
	  } else if (fabs(ybincenter-1.85) < 0.01 ) {
	    theBin = histJpsi16_21.FindBin(ybincenter,ptbincenter);
	    histJpsi16_21.SetBinContent(theBin,fitted);
	    histJpsi16_21.SetBinError(theBin,error);  
	    histRelVar16_21.SetBinContent(theBin,relvar);
	  } else {
	    theBin = histJpsi21_24.FindBin(ybincenter,ptbincenter);
	    histJpsi21_24.SetBinContent(theBin,fitted);
	    histJpsi21_24.SetBinError(theBin,error); 
	    histRelVar21_24.SetBinContent(theBin,relvar);
	    // histErr.SetBinContent(theBin,error);
	    // histPul.SetBinContent(theBin,(fitted-trueMC)/error);
	    // chi2i += pow((fitted-trueMC)/error,2);
	  }
	  i++;
	}
      }  

      theFile.close();
      theFile2.close();
      if (maxBetween >= 2) theFile3.close();
      if (maxBetween == 3) theFile4.close();
    }
  }

  lista.close();
  lista2.close();
  if (maxBetween >= 2) lista3.close();
  if (maxBetween == 3) lista4.close();

  cutstring = "results/" + fileNameBase + whichTypeJpsi + "_summary.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histJpsi00_09.Write();
  histJpsi09_12.Write();
  histJpsi12_16.Write();
  histJpsi16_21.Write();
  histJpsi21_24.Write();
  histRelVar00_09.Write();
  histRelVar09_12.Write();
  histRelVar12_16.Write();
  histRelVar16_21.Write();
  histRelVar21_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
