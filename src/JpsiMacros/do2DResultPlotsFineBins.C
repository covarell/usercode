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

  static const unsigned int nbinspt1 = 12;
  static const unsigned int nbinspt2 = 8;
  static const unsigned int nbinspt3 = 13;
  static const unsigned int nbinspt4 = 22;
  static const unsigned int nbinspt5 = 9;

  double pt1binlimits[nbinspt1+1] = {6.5,7.,8.,9.,10.,11.,12.,13.5,15.,18.,30.,45.,70.};
  double pt2binlimits[nbinspt2+1] = {6.0,6.5,8.,9.,10.,12.,15.,30.,45.};
  double pt3binlimits[nbinspt3+1] = {5.0,6.0,6.5,7.,7.5,8.,8.5,9.,10.,11.,12.,15.,30.,45.};
  double pt4binlimits[nbinspt4+1] = {4.,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.,6.25,6.5,6.75,7.,7.25,7.5,8.,8.5,9.,10.,11.,12.,15.,30.};
  double pt5binlimits[nbinspt5+1] = {4.,4.5,5.5,6.5,8.,9.,10.,12.,15.,30.};
  
  TH2F histBfracJpsi00_09("histBfracJpsi00_09","histBfracJpsi",5,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histBfracJpsi09_12("histBfracJpsi09_12","histBfracJpsi",5,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histBfracJpsi12_16("histBfracJpsi12_16","histBfracJpsi",5,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histBfracJpsi16_21("histBfracJpsi16_21","histBfracJpsi",5,ybinlimits,nbinspt4,pt4binlimits);
  TH2F histBfracJpsi21_24("histBfracJpsi21_24","histBfracJpsi",5,ybinlimits,nbinspt5,pt5binlimits);

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
	    theBin = histBfracJpsi00_09.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi00_09.SetBinContent(theBin,fitted);
	    histBfracJpsi00_09.SetBinError(theBin,error);
	  } else if (fabs(ybincenter-1.05) < 0.01 ) {
	    theBin = histBfracJpsi09_12.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi09_12.SetBinContent(theBin,fitted);
	    histBfracJpsi09_12.SetBinError(theBin,error);   
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histBfracJpsi12_16.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi12_16.SetBinContent(theBin,fitted);
	    histBfracJpsi12_16.SetBinError(theBin,error); 
	  } else if (fabs(ybincenter-1.85) < 0.01 ) {
	    theBin = histBfracJpsi16_21.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi16_21.SetBinContent(theBin,fitted);
	    histBfracJpsi16_21.SetBinError(theBin,error);  
	  } else {
	    theBin = histBfracJpsi21_24.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi21_24.SetBinContent(theBin,fitted);
	    histBfracJpsi21_24.SetBinError(theBin,error); 
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

  cutstring = "pictures/" + fileNameBase + whichTypeJpsi + "_results.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histBfracJpsi00_09.Write();
  histBfracJpsi09_12.Write();
  histBfracJpsi12_16.Write();
  histBfracJpsi16_21.Write();
  histBfracJpsi21_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
