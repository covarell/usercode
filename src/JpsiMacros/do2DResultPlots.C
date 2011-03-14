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

void do2DResultPlots(string whichTypeJpsi = "BP", string whichTypePsip = "BJ", string fileNameBase = "results") {

  double ybinlimits[4] = {0.0,1.2,1.6,2.4};

  static const unsigned int nbinspt1 = 6;
  static const unsigned int nbinspt2 = 4;
  static const unsigned int nbinspt3 = 6;

  double pt1binlimits[nbinspt1+1] = {6.5,8.,9.,10.,12.,15.,30.};
  double pt2binlimits[nbinspt2+1] = {5.5,6.5,8.,10.,30.};
  double pt3binlimits[nbinspt3+1] = {4.5,5.5,6.5,8.,10.,12.,30.};
  
  TH2F histBfracJpsi00_12("histBfracJpsi00_12","histBfracJpsi",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histBfracJpsi12_16("histBfracJpsi12_16","histBfracJpsi",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histBfracJpsi16_24("histBfracJpsi16_24","histBfracJpsi",3,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histBfracPsip00_12("histBfracPsip00_12","histBfracPsip",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histBfracPsip12_16("histBfracPsip12_16","histBfracPsip",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histBfracPsip16_24("histBfracPsip16_24","histBfracPsip",3,ybinlimits,nbinspt3,pt3binlimits);
  // TH2F histErr("histErr","histErr",2,ybinlimits,nbinspt,ptbinlimits);
  // TH2F histPul("histPul","histPul",2,ybinlimits,nbinspt,ptbinlimits);

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
	  if (fabs(ybincenter-0.6) < 0.01 ) {
	    theBin = histBfracJpsi00_12.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi00_12.SetBinContent(theBin,fitted);
	    histBfracJpsi00_12.SetBinError(theBin,error);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histBfracJpsi12_16.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi12_16.SetBinContent(theBin,fitted);
	    histBfracJpsi12_16.SetBinError(theBin,error);  
	  } else {
	    theBin = histBfracJpsi16_24.FindBin(ybincenter,ptbincenter);
	    histBfracJpsi16_24.SetBinContent(theBin,fitted);
	    histBfracJpsi16_24.SetBinError(theBin,error); 
	    // histErr.SetBinContent(theBin,error);
	    // histPul.SetBinContent(theBin,(fitted-trueMC)/error);
	    // chi2i += pow((fitted-trueMC)/error,2);
	  }
	  i++;
	  
	} else if (!strcmp(theType,whichTypePsip.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin + 0.01)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          cout << "TEST : " << ptbincenter << " " << ybincenter << " " << theType << " " << trueMC << " " << fitted << " " << error << endl;
	  if (fabs(ybincenter-0.6) < 0.01 ) {
	    theBin = histBfracPsip00_12.FindBin(ybincenter,ptbincenter);
	    histBfracPsip00_12.SetBinContent(theBin,fitted);
	    histBfracPsip00_12.SetBinError(theBin,error);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histBfracPsip12_16.FindBin(ybincenter,ptbincenter);
	    histBfracPsip12_16.SetBinContent(theBin,fitted);
	    histBfracPsip12_16.SetBinError(theBin,error);  
	  } else {
	    theBin = histBfracPsip16_24.FindBin(ybincenter,ptbincenter);
	    histBfracPsip16_24.SetBinContent(theBin,fitted);
	    histBfracPsip16_24.SetBinError(theBin,error); 
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

  cutstring = "pictures/" + fileNameBase + whichTypeJpsi + whichTypePsip + "_results.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histBfracJpsi00_12.Write();
  histBfracJpsi12_16.Write();
  histBfracJpsi16_24.Write();
  histBfracPsip00_12.Write();
  histBfracPsip12_16.Write();
  histBfracPsip16_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
