// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2F.h>

void produceOtherMikeTables(int LHCsqrts = 7, string sample = "gg") {

  cout << "FIRST RUN source changeCinZero.csh!" << endl;

  string Sources[6] = {"","","","","","",""};
  if (sample == "gg") {
    Sources[0] = "Resummation";
    Sources[1] = "TopMass";
  }
  if (sample == "vbf") {
    Sources[0] = "PDF-VBF";
    Sources[1] = "scale-VBF";
  }
  if (sample == "zz") {
    Sources[0] = "SingleZ";
    Sources[1] = "PDF";
    Sources[2] = "scale";
  }
  // if (sample == "zz" && LHCsqrts == 8) {
  //  Sources[0] = "UnbRegion";
  // }

  ofstream *mikeFile;
  char fileName[200];

  char thePar[10];
  char equalS[1];
  char dashS[1];
  char pmS[3];
  char theLimit1[10];
  char theLimit2[10];
  float fitted, error;

  const int nPars = 7;
  string parNames[nPars] = {"m","n","n2","bb","bb2","fexp","T"};
  float parVal[nPars];
  float upVar[nPars];
  float downVar[nPars];

  // Inizia
  // sprintf(fileName,"MELA/scripts/chris/paramsCJLST_%s_%dTeV_Default.txt",sample.c_str(),LHCsqrts);
  sprintf(fileName,"text/paramsPTOverMCJLST_%s_%dTeV_Default.txt",sample.c_str(),LHCsqrts);
  ifstream theFile(fileName);

  sprintf(fileName,"text/paramShifts_%s_%dTeV_Statistics.txt",sample.c_str(),LHCsqrts);
  mikeFile = new ofstream(fileName);
  
  while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
    cout << thePar << " " << fitted << " " << error << endl;
    for (int ii = 0; ii < nPars; ii++) {
      string theString = parNames[ii]+"up";
      // string theString = parNames[ii];
      if (!strcmp(thePar,theString.c_str())) {
	parVal[ii] = fitted; 
	upVar[ii] = fitted+error;
	downVar[ii] = fabs(fitted-error);
      }
    }
  }

  theFile.close();

  for (int ii = 0; ii < nPars; ii++) {
    *mikeFile << parNames[ii].c_str() << " " << parVal[ii] << " " << upVar[ii] << " " << downVar[ii] << endl; 
  }
   

  for (int jj = 0; jj < 6; jj++) {
    if (Sources[jj] != "") {
      
      cout << "Systematics n. " << jj+1 << " : " << Sources[jj] << endl;

      // sprintf(fileName,"MELA/scripts/chris/paramsCJLST_%s_%dTeV_%s.txt",sample.c_str(),LHCsqrts,Sources[jj].c_str());
      sprintf(fileName,"text/paramsPTOverMCJLST_%s_%dTeV_%s.txt",sample.c_str(),LHCsqrts,Sources[jj].c_str());
      ifstream theFile(fileName);

      sprintf(fileName,"text/paramShifts_%s_%dTeV_%s.txt",sample.c_str(),LHCsqrts,Sources[jj].c_str());
      mikeFile = new ofstream(fileName);

      while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	cout << thePar << " " << fitted << " " << error << endl;
	for (int ii = 0; ii < nPars; ii++) {
	  if (!strcmp(thePar,parNames[ii].c_str())) {
	    upVar[ii] = fitted;
	    downVar[ii] = fabs(2*parVal[ii] - fitted);
	  }
	}
      }

      for (int ii = 0; ii < nPars; ii++) {
	*mikeFile << parNames[ii].c_str() << " " << parVal[ii] << " " << upVar[ii] << " " << downVar[ii] << endl; 
      }
   
      theFile.close();

    }
  }


  return;
}
