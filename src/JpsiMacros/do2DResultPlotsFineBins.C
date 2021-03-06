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

void do2DResultPlotsFineBins(string whichTypeJpsi = "BP", string fileNameBase = "results") {

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

  gROOT->ProcessLine(".! rm -f lista");
  char theCommand[200];
  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  
  ifstream lista("lista");
  char fileName[300];  
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  float trueMC, fitted, error; 
  int i = 0;       
  string cutstring;
  float chi2i = 0.;  
  float theMaximumPt = -1.;
  float theMinimumPt = 9000000.;
  int theBin = 0;
  // Inizia

  while (lista >> fileName) {

    cutstring = "results/" + fileNameBase + "_pT%f-%f_y%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
	
      while (theFile >> theType >> trueMC >> fitted >> error) {
	
	if (!strcmp(theType,whichTypeJpsi.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin + 0.01)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          cout << "TEST : " << ptbincenter << " " << ybincenter << " " << theType << " " << trueMC << " " << fitted << " " << error << endl;
	  if (fabs(ybincenter-0.45) < 0.01 ) {
	    theBin = histJpsi00_09.FindBin(ybincenter,ptbincenter);
	    histJpsi00_09.SetBinContent(theBin,fitted);
	    histJpsi00_09.SetBinError(theBin,error);
	  } else if (fabs(ybincenter-1.05) < 0.01 ) {
	    theBin = histJpsi09_12.FindBin(ybincenter,ptbincenter);
	    histJpsi09_12.SetBinContent(theBin,fitted);
	    histJpsi09_12.SetBinError(theBin,error);   
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histJpsi12_16.FindBin(ybincenter,ptbincenter);
	    histJpsi12_16.SetBinContent(theBin,fitted);
	    histJpsi12_16.SetBinError(theBin,error); 
	  } else if (fabs(ybincenter-1.85) < 0.01 ) {
	    theBin = histJpsi16_21.FindBin(ybincenter,ptbincenter);
	    histJpsi16_21.SetBinContent(theBin,fitted);
	    histJpsi16_21.SetBinError(theBin,error);  
	  } else {
	    theBin = histJpsi21_24.FindBin(ybincenter,ptbincenter);
	    histJpsi21_24.SetBinContent(theBin,fitted);
	    histJpsi21_24.SetBinError(theBin,error); 
	    // histErr.SetBinContent(theBin,error);
	    // histPul.SetBinContent(theBin,(fitted-trueMC)/error);
	    // chi2i += pow((fitted-trueMC)/error,2);
	  }
	  i++;
	}
      }  

      theFile.close();
    }
  }

  lista.close();

  cutstring = "results/" + fileNameBase + whichTypeJpsi + "_summary.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histJpsi00_09.Write();
  histJpsi09_12.Write();
  histJpsi12_16.Write();
  histJpsi16_21.Write();
  histJpsi21_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
