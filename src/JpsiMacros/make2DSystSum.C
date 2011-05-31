// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

void make2DSystSum() {

  string theDir = "/afs/cern.ch/user/c/covarell/public/html/quarkonia/JpsiPsipBfrac/";
  string filenames[5] = {"systLifetimeModel/results2DALLBJ_BP_RATIO_summary.root",
			 "systPrimaryV/results2DALLBJ_BP_RATIO_summary.root",
			 "systResoModel/results2DALLBJ_BP_RATIO_summary.root",
			 "systAlignment/results2DALLBJ_BP_RATIO_summary.root",
			 "systSidebands/results2DALLBJ_BP_RATIO_summary.root"};
  
  // string histonames[5] = {"histRelVar00_09",
  //			  "histRelVar09_12",
  //			  "histRelVar12_16",
  //     		  "histRelVar16_21",
  //		          "histRelVar21_24"};

  string histonames[3] = {"histRelVarRatio00_12",
			  "histRelVarRatio12_16",
			  "histRelVarRatio16_24"};
  
 
  TFile* f1;
  TH2D* h2;
  TH2D* newhistos[5];
  string totname = theDir + filenames[0];
  f1 = TFile::Open(totname.c_str());
  for (unsigned int jj = 0; jj < 3; jj++) {
    newhistos[jj] = (TH2D*)(f1->Get(histonames[jj].c_str()))->Clone();
    newhistos[jj]->SetName(histonames[jj].c_str());
  }

  for (unsigned int i = 1; i < 5; i++) {
    totname = theDir + filenames[i];
    f1 = TFile::Open(totname.c_str());
    for (unsigned int jj = 0; jj < 3; jj++) {
      h2 = (TH2D*)f1->Get(histonames[jj].c_str());
      for (unsigned int j = 1; j <= h2->GetNbinsX(); j++) {
	for (unsigned int k = 1; k <= h2->GetNbinsY(); k++) {
	  if (newhistos[jj]->GetBinContent(j,k) > 0.0001) {
	    float newBinCont;
	    if (i == 1) newBinCont = pow(newhistos[jj]->GetBinContent(j,k),2) + pow(h2->GetBinContent(j,k),2);
	    else newBinCont = newhistos[jj]->GetBinContent(j,k) + pow(h2->GetBinContent(j,k),2);
	    cout << i << " " << jj << " " << j << " " << k << " " << newhistos[jj]->GetBinContent(j,k) << " " << h2->GetBinContent(j,k) << endl;
	    newhistos[jj]->SetBinContent(j,k,newBinCont);
	  }
	}
      }
    }
  }

  for (unsigned int jj = 0; jj < 3; jj++) {
    for (unsigned int j = 1; j <= newhistos[jj]->GetNbinsX(); j++) {
      for (unsigned int k = 1; k <= newhistos[jj]->GetNbinsY(); k++) {
	if (newhistos[jj]->GetBinContent(j,k) > 0.0001) {
	  float newBinCont = sqrt(newhistos[jj]->GetBinContent(j,k));
	  newhistos[jj]->SetBinContent(j,k,newBinCont);
	}
      }
    }
  }
  
  TFile out("results/results2DALLRATIO_summary.root","RECREATE");
  for (unsigned int jj = 0; jj < 3; jj++) {
    newhistos[jj]->Write();
  }
  out.Close();
  return; 
}
