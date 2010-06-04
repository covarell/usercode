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

void doMeanPlots(string fileNameBase = "results") {

  static const unsigned int nbinspt = 6;

  double ptbincenters[nbinspt+1] = {3.0,6.0,8.0,9.0,10.0,13.0,30.0};
  double ybincenters[3] = {0.0,1.4,2.5};
  
  TH2F meanPtGG("meanPtGG","meanPtGG",2,ybincenters,nbinspt,ptbincenters);
  TH2F meanPt2GG("meanPt2GG","meanPt2GG",2,ybincenters,nbinspt,ptbincenters);
  TH2F meanPtGT("meanPtGT","meanPtGT",2,ybincenters,nbinspt,ptbincenters);
  TH2F meanPt2GT("meanPt2GT","meanPt2GT",2,ybincenters,nbinspt,ptbincenters);

  gROOT->ProcessLine(".! rm -f lista");
  char theCommand[50];
  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  
  ifstream lista("lista");
  char fileName[200];  
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  float meanP, meanP2, dummy; 
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
	
      while (theFile >> theType >> meanP >> meanP2 >> dummy) {
	
	if (!strcmp(theType,"GG")) {
	  
	  float ybincenter = (etamax + etamin)/2. ; 
	  int theBin = meanPtGG.FindBin(ybincenter,meanP);
	  meanPtGG.SetBinContent(theBin,meanP);
	  meanPt2GG.SetBinContent(theBin,meanP2);
	  
	}
      }  

      theFile.close();
    }
  }

  lista.close();
  ifstream lista2("lista");
  while (lista2 >> fileName) {

    cutstring = "results/" + fileNameBase + "_pT%f-%f_eta%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
	
      while (theFile >> theType >> meanP >> meanP2 >> dummy) {
	
	if (!strcmp(theType,"GT")) {
	  
	  float ybincenter = (etamax + etamin)/2. ; 
	  int theBin = meanPtGT.FindBin(ybincenter,meanP);
	  meanPtGT.SetBinContent(theBin,meanP);
	  meanPt2GT.SetBinContent(theBin,meanP2);
	  
	}
      }  

      theFile.close();
    }
  }


  cutstring = "pictures/" + fileNameBase + "_means.root"; 
  TFile f55(cutstring.c_str(),"RECREATE");  
  meanPtGG.Write();
  meanPt2GG.Write();
  meanPtGT.Write();
  meanPt2GT.Write();
  f55.Close();   
  cout << " Done " << endl;

  return;
}
