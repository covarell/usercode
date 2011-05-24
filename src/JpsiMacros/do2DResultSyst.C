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

void do2DResultSyst(string whichTypeJpsi = "BP", string whichTypePsip = "BJ", string fileNameBase = "results", string fileNameComp = "results2", int maxBetween = 1, string fileNameComp2 = "none", string fileNameComp3 = "none") {

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
  TH2F histRelVarPsip00_12("histRelVarPsip00_12","histPsip - relative variation",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histRelVarPsip12_16("histRelVarPsip12_16","histPsip - relative variation",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histRelVarPsip16_24("histRelVarPsip16_24","histPsip - relative variation",3,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histRatio00_12("histRatio00_12","histRatio",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histRatio12_16("histRatio12_16","histRatio",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histRatio16_24("histRatio16_24","histRatio",3,ybinlimits,nbinspt3,pt3binlimits);
  TH2F histRelVarRatio00_12("histRelVarRatio00_12","histRatio - relative variation",3,ybinlimits,nbinspt1,pt1binlimits);
  TH2F histRelVarRatio12_16("histRelVarRatio12_16","histRatio - relative variation",3,ybinlimits,nbinspt2,pt2binlimits);
  TH2F histRelVarRatio16_24("histRelVarRatio16_24","histRatio - relative variation",3,ybinlimits,nbinspt3,pt3binlimits);
  // TH2F histErr("histErr","histErr",2,ybinlimits,nbinspt,ptbinlimits);
  // TH2F histPul("histPul","histPul",2,ybinlimits,nbinspt,ptbinlimits);

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
  float trueMC, fitted, error, fittedJ, errorJ; 
  float fitted2, error2, fittedJ2, errorJ2;
  float fitted3, fittedJ3; 
  float fitted4, fittedJ4; 
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
	if (maxBetween >= 2) theFile3 >> theType2 >> trueMC >> fitted3 >> error2;
	if (maxBetween == 3) theFile4 >> theType2 >> trueMC >> fitted4 >> error2;

	if (!strcmp(theType,whichTypeJpsi.c_str())) {
	  
	  float ptbincenter = (ptmax + ptmin + 0.01)/2. ;
	  float ybincenter = (etamax + etamin)/2. ; 
          cout << "TEST : " << ptbincenter << " " << ybincenter << " " << theType << " " << trueMC << " " << fitted << " " << error << endl;
	  fittedJ = fitted;   errorJ = error;
	  fittedJ2 = fitted2;   errorJ2 = error2;
	  if (maxBetween >= 2) fittedJ3 = fitted3;
	  if (maxBetween == 3) fittedJ4 = fitted4;
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
	  float relvarP = 100.0*fabs(fitted-fitted2)/fitted2;
	  if (maxBetween >= 2) {
	    float relvar2P = 100.0*fabs(fitted2-fitted3)/fitted2;
	    if (relvar2P > relvarP) {
	      relvarP = relvar2P;   fitted = fitted3;
	    }
	  }
	  if (maxBetween == 3) {
	    float relvar3P = 100.0*fabs(fitted2-fitted4)/fitted2;
	    if (relvar3P > relvarP) {
	      relvarP = relvar3P;   fitted = fitted4;
	    }
	  }
	  float ratio = fitted/fittedJ;
	  float errratio = ratio*sqrt(pow(error/fitted,2) + pow(errorJ/fittedJ,2));
	  float ratio2 = fitted2/fittedJ2;
	  float relvarR = 100.0*fabs(ratio-ratio2)/ratio2;
	  if (maxBetween >= 2) {
	    float ratio3 = fitted3/fittedJ3;
	    float relvar2R = 100.0*fabs(ratio2-ratio3)/ratio2;
	    if (relvar2R > relvarR) {
	      relvarR = relvar2R;    ratio = ratio3;
	    }
	  }
	  if (maxBetween == 3) {
	    float ratio4 = fitted4/fittedJ4;
	    float relvar3R = 100.0*fabs(ratio2-ratio4)/ratio2;
	    if (relvar3R > relvarR) {
	      relvarR = relvar3R;    ratio = ratio4;
	    }
	  }
	  if (fabs(ybincenter-0.6) < 0.01 ) {
	    theBin = histPsip00_12.FindBin(ybincenter,ptbincenter);
	    histPsip00_12.SetBinContent(theBin,fitted);
	    histPsip00_12.SetBinError(theBin,error);
	    histRatio00_12.SetBinContent(theBin,ratio);
	    histRatio00_12.SetBinError(theBin,errratio);
	    histRelVarPsip00_12.SetBinContent(theBin,relvarP);
	    histRelVarRatio00_12.SetBinContent(theBin,relvarR);
	  } else if (fabs(ybincenter-1.4) < 0.01 ) {
	    theBin = histPsip12_16.FindBin(ybincenter,ptbincenter);
	    histPsip12_16.SetBinContent(theBin,fitted);
	    histPsip12_16.SetBinError(theBin,error);
	    histRatio12_16.SetBinContent(theBin,ratio);
	    histRatio12_16.SetBinError(theBin,errratio);
	    histRelVarPsip12_16.SetBinContent(theBin,relvarP);
	    histRelVarRatio12_16.SetBinContent(theBin,relvarR);
	  } else {
	    theBin = histPsip16_24.FindBin(ybincenter,ptbincenter);
	    histPsip16_24.SetBinContent(theBin,fitted);
	    histPsip16_24.SetBinError(theBin,error); 
	    histRatio16_24.SetBinContent(theBin,ratio);
	    histRatio16_24.SetBinError(theBin,errratio);
	    histRelVarPsip16_24.SetBinContent(theBin,relvarP);
	    histRelVarRatio16_24.SetBinContent(theBin,relvarR);
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
  histRelVarPsip00_12.Write();
  histRelVarPsip12_16.Write();
  histRelVarPsip16_24.Write();
  histRelVarRatio00_12.Write();
  histRelVarRatio12_16.Write();
  histRelVarRatio16_24.Write();
  f55.Close();   
  cout << "chi2 = " << chi2i/6. << endl;

  return;
}
