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

  static const unsigned int nbinspt1 = 9;
  static const unsigned int nbinspt2 = 7;
  static const unsigned int nbinspt3 = 7;

  double pt1binlimits[nbinspt1+1] = {6.5,8.,9.,10.,11.,12.,13.5,15.,18.,30.};
  double pt2binlimits[nbinspt2+1] = {5.5,6.5,8.,9.,10.,12.,15.,30.};
  double pt3binlimits[nbinspt3+1] = {5.5,6.5,8.,9.,10.,12.,15.,30.};
  
  TH2F histJpsi00_12("histJpsi00_12","histJpsi",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histJpsi12_16("histJpsi12_16","histJpsi",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histJpsi16_24("histJpsi16_24","histJpsi",3,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histPsip00_12("histPsip00_12","histPsip",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histPsip12_16("histPsip12_16","histPsip",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histPsip16_24("histPsip16_24","histPsip",3,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histRatio00_12("histRatio00_12","histRatio",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histRatio12_16("histRatio12_16","histRatio",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histRatio16_24("histRatio16_24","histRatio",3,ybinlimits,nbinspt3,pt3binlimits);
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
  float trueMC, fitted, error, fittedJ, errorJ; 
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
	  fittedJ = fitted;   errorJ = error;
	  if (fabs(ybincenter-0.6) < 0.01 ) {
	    theBin = histJpsi00_12.FindBin(ybincenter,ptbincenter);
	    histJpsi00_12.SetBinContent(theBin,fitted);
	    histJpsi00_12.SetBinError(theBin,error);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histJpsi12_16.FindBin(ybincenter,ptbincenter);
	    histJpsi12_16.SetBinContent(theBin,fitted);
	    histJpsi12_16.SetBinError(theBin,error);  
	  } else {
	    theBin = histJpsi16_24.FindBin(ybincenter,ptbincenter);
	    histJpsi16_24.SetBinContent(theBin,fitted);
	    histJpsi16_24.SetBinError(theBin,error); 
	    // histErr.SetBinContent(theBin,error);
	    // histPul.SetBinContent(theBin,(fitted-trueMC)/error);
	    // chi2i += pow((fitted-trueMC)/error,2);
	  }
	  i++;
	  
	} else if (!strcmp(theType,whichTypePsip.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin + 0.01)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          cout << "TEST : " << ptbincenter << " " << ybincenter << " " << theType << " " << trueMC << " " << fitted << " " << error << endl;
	  float ratio = fitted/fittedJ;
	  float errratio = ratio*sqrt(pow(error/fitted,2) + pow(errorJ/fittedJ,2));
	  if (fabs(ybincenter-0.6) < 0.01 ) {
	    theBin = histPsip00_12.FindBin(ybincenter,ptbincenter);
	    histPsip00_12.SetBinContent(theBin,fitted);
	    histPsip00_12.SetBinError(theBin,error);
	    histRatio00_12.SetBinContent(theBin,ratio);
	    histRatio00_12.SetBinError(theBin,errratio);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histPsip12_16.FindBin(ybincenter,ptbincenter);
	    histPsip12_16.SetBinContent(theBin,fitted);
	    histPsip12_16.SetBinError(theBin,error);
	    histRatio12_16.SetBinContent(theBin,ratio);
	    histRatio12_16.SetBinError(theBin,errratio);
	  } else {
	    theBin = histPsip16_24.FindBin(ybincenter,ptbincenter);
	    histPsip16_24.SetBinContent(theBin,fitted);
	    histPsip16_24.SetBinError(theBin,error); 
	    histRatio16_24.SetBinContent(theBin,ratio);
	    histRatio16_24.SetBinError(theBin,errratio);
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

  cutstring = "results/" + fileNameBase + whichTypeJpsi + "_" + whichTypePsip + "_RATIO_summary.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histJpsi00_12.Write();
  histJpsi12_16.Write();
  histJpsi16_24.Write();
  histPsip00_12.Write();
  histPsip12_16.Write();
  histPsip16_24.Write();
  histRatio00_12.Write();
  histRatio12_16.Write();
  histRatio16_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
