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

void do2DResultPlots(string whichType = "GG", string fileNameBase = "results") {

  static const unsigned int nbinspt = 6;

  double ptbincenters[nbinspt+1] = {3.0,6.0,8.0,9.0,10.0,13.0,30.0};
  double ybincenters[3] = {0.0,1.4,2.5};
  
  TH2F histRes("histRes","histRes",2,ybincenters,nbinspt,ptbincenters);
  TH2F histErr("histErr","histErr",2,ybincenters,nbinspt,ptbincenters);
  TH2F histPul("histPul","histPul",2,ybincenters,nbinspt,ptbincenters);

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

  // Inizia

  while (lista >> fileName) {

    cutstring = "results/" + fileNameBase + "_pT%f-%f_eta%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
	
      while (theFile >> theType >> trueMC >> fitted >> error) {
	
	if (!strcmp(theType,whichType.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          cout << "TEST : " << ptbincenter << " " << ybincenter << " " << theType << " " << trueMC << " " << fitted << " " << error << endl;
	  int theBin = histRes.FindBin(ybincenter,ptbincenter);
	  histRes.SetBinContent(theBin,fitted);
          histRes.SetBinError(theBin,error);
	  histErr.SetBinContent(theBin,error);
	  histPul.SetBinContent(theBin,(fitted-trueMC)/error);
	  chi2i += pow((fitted-trueMC)/error,2);
	  i++;
	  
	}
      }  

      theFile.close();
    }
  }

  lista.close();

  cutstring = "pictures/" + fileNameBase + whichType + "_results.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  histRes.Write();
  histErr.Write();
  histPul.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/12. << endl;

  return;
}
