// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

void doDeviationPlot(string whichType = "GG", string fileNameBase = "results", float minMax = 500.) {

  static const unsigned int nbinspt = 7;

  double ptbincenters[nbinspt] = {nbinspt*20.};
  double ptbinerrors[nbinspt] = {nbinspt*0.};
  double ycenters1[nbinspt] = {nbinspt*-999.};
  double yerrors1[nbinspt] = {nbinspt*0.};
  double ycenters2[nbinspt] = {nbinspt*-999.};
  double yerrors2[nbinspt] = {nbinspt*0.};
  
  gROOT->ProcessLine(".! rm -f lista");
  char theCommand[50];
  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  
  ifstream lista("lista");
  char fileName[200];  
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  float trueMC, fitted, error; 
  int i = 0; 
  int j = 0;
  string cutstring;

  // Inizia

  while (!lista.eof()) {

    lista >> fileName;
    cutstring = "results/" + fileNameBase + "_pT%f-%f_eta%f-%f.txt";
    cout << "1 " << fileName << endl;  
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
      cout << "2 " << fileName << endl;  
      if (etamin == 0.0) {
	
	while (theFile >> theType >> trueMC >> fitted >> error) {
	
          cout << i << " " << theType << " " << trueMC << endl;
	  if (!strcmp(theType,whichType.c_str())) {

	    ptbincenters[i] = (ptmax + ptmin)/2. ;
	    ptbinerrors[i] = (ptmax - ptmin)/2. ; 
	    ycenters1[i] = fitted - trueMC;
	    yerrors1[i] = error;
	    i++;
	    cout << i << " " << fileName << " " << etamin << " " << etamax << endl;
	  }
	}  
      } else if (etamax == 2.5) {
	while (theFile >> theType >> trueMC >> fitted >> error) {

	  if (!strcmp(theType,whichType.c_str())) {
	    ycenters2[j] = fitted - trueMC;
	    yerrors2[j] = error;
	    j++;
	  }
	} 
      }

      theFile.close();
    }
  }

  lista.close();

  for (int aa = 0; aa < nbinspt; aa++) {
    cout << ptbincenters[aa] << " " << ptbinerrors[aa] << " " << ycenters1[aa] << " " << yerrors1[aa] << " " << ycenters2[aa] << " " << yerrors2[aa] << " " << endl;
  }

  TGraphErrors *gpull1 = new TGraphErrors(nbinspt,ptbincenters,ycenters1,ptbinerrors,yerrors1);

  TCanvas c3;
  c3.cd();
  gpull1->SetTitle("Pull of fit results - barrel");
  gpull1->SetMarkerStyle(20);
  gpull1->SetMarkerColor(kRed);
  gpull1->SetMaximum(minMax);
  gpull1->SetMinimum(-minMax);
  gpull1->Draw("AP");

  cutstring = whichType + "pull_barrel.gif";
  c3.SaveAs(cutstring.c_str());

  TGraphErrors *gpull2 = new TGraphErrors(nbinspt,ptbincenters,ycenters2,ptbinerrors,yerrors2);

  TCanvas c4;
  c4.cd();
  gpull2->SetTitle("Pull of fit results GG - endcap");
  gpull2->SetMarkerStyle(20);
  gpull2->SetMarkerColor(kBlue);
  gpull2->SetMaximum(minMax);
  gpull2->SetMinimum(-minMax);
  gpull2->Draw("AP");

  cutstring = whichType + "pull_endcap.gif";
  c4.SaveAs(cutstring.c_str());
        
  return;
}
