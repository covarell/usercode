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

void produceMikeTables(int LHCsqrts = 7, string sample = "gg") {

  cout << "FIRST RUN source changeCinZero.csh!" << endl;

  string Sources[6] = {"","","","","","",""};
  if (sample == "gg") {
    Sources[0] = "Resummation";
    Sources[1] = "TopMass";
  }
  if (sample == "zz" && LHCsqrts == 7) {
    Sources[0] = "SingleZ";
  }
  if (sample == "zz" && LHCsqrts == 8) {
    Sources[0] = "UnbRegion";
  }

  ofstream *mikeFile;
  char fileName[200];

  sprintf(fileName,"text/paramShifts_%s_%dTeV.txt",sample.c_str(),LHCsqrts);
  mikeFile = new ofstream(fileName);

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
  sprintf(fileName,"text/paramsCJLST_%s_%dTeV_default.txt",sample.c_str(),LHCsqrts);
  ifstream theFile(fileName);
  
  while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
    cout << thePar << " " << fitted << " " << error << endl;
    for (int ii = 0; ii < nPars; ii++) {
      if (!strcmp(thePar,parNames[ii].c_str())) {
	parVal[ii] = fitted; 
	upVar[ii] = error*error;
	downVar[ii] = error*error;
      }
    }
  }
   
  theFile.close();

  for (int jj = 0; jj < 6; jj++) {
    if (Sources[jj] != "") {
      
      cout << "Systematics n. " << jj+1 << " : " << Sources[jj] << endl;

      sprintf(fileName,"text/paramsCJLST_%s_%dTeV_%s.txt",sample.c_str(),LHCsqrts,Sources[jj].c_str());
      ifstream theFile(fileName);
      while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	cout << thePar << " " << fitted << " " << error << endl;
	for (int ii = 0; ii < nPars; ii++) {
	  if (!strcmp(thePar,parNames[ii].c_str())) {
	    if (fitted > parVal[ii]) upVar[ii] += (fitted - parVal[ii])*(fitted - parVal[ii]);
	    if (fitted < parVal[ii]) downVar[ii] += (fitted - parVal[ii])*(fitted - parVal[ii]);
	  }
	}
      }

      theFile.close();
    }
  }

  for (int ii = 0; ii < nPars; ii++) {
    upVar[ii] = sqrt(upVar[ii]);
    downVar[ii] = sqrt(downVar[ii]);
    *mikeFile << parNames[ii].c_str() << " " << parVal[ii] << " " << parVal[ii]+upVar[ii] << " " << parVal[ii]-downVar[ii] << endl; 
  }

  return;
}
