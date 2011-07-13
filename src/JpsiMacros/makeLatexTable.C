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

void makeLatexTable(string rootfile = 'myfile.root') {
// string systfile = 'mysyst.root', string systfile2 = 'mysyst2.root') {

 //  string filenames[5] = {"~/public/html/quarkonia/Aeff2/Aeff_New_y0_sg_jpsi.root",
// 			 "~/public/html/quarkonia/Aeff2/Aeff_New_y1_sg_jpsi.root",
// 			 "~/public/html/quarkonia/Aeff2/Aeff_New_y2_sg_jpsi.root",
// 			 "~/public/html/quarkonia/Aeff2/Aeff_New_y3_sg_jpsi.root",
// 			 "~/public/html/quarkonia/Aeff2/Aeff_New_y4_sg_jpsi.root"};

  string filenames[3] = {"~/public/html/quarkonia/Aeff2/Aeff_New_y2s0_sg_psi2s.root",
			 "~/public/html/quarkonia/Aeff2/Aeff_New_y2s1_sg_psi2s.root",
			 "~/public/html/quarkonia/Aeff2/Aeff_New_y2s2_sg_psi2s.root"};

  string histonames[5] = {"histJpsi00_09",
			  "histJpsi09_12",
			  "histJpsi12_16",
			  "histJpsi16_21",
			  "histJpsi21_24"};
  
 //  string histosyst[5]  = {"histRelVar00_09",
// 			  "histRelVar09_12",
// 			  "histRelVar12_16",
// 			  "histRelVar16_21",
// 			  "histRelVar21_24"};


//   string histonames[3] = {"histPsip00_12",
// 			  "histPsip12_16",
// 			  "histPsip16_24"};

//   string histonames2[3] = {"histRatio00_12",
// 			  "histRatio12_16",
// 			  "histRatio16_24"};
  
//   string histosyst[3] = {"histRelVarPsip00_12",
// 			  "histRelVarPsip12_16",
// 			  "histRelVarPsip16_24"};

//   string histosyst2[3] = {"histRelVarRatio00_12",
// 			  "histRelVarRatio12_16",
// 			  "histRelVarRatio16_24"};


  // TFile f(rootfile.c_str());
  // TFile fs(systfile.c_str());
  // TFile fs2(systfile2.c_str());

  TH2D* histo, histo2, histo3, histo4;
	
  for (unsigned int i = 0; i < 3; i++) {
    TFile f(filenames[i].c_str());
    // histo = (TH2D*)f.Get(histonames[i].c_str());
    histo = (TH2D*)f.Get("Aeff_fitted");
    // histo2 = (TH2D*)fs.Get(histosyst[i].c_str());
    // histo3 = (TH2D*)f.Get(histonames2[i].c_str());
    // histo4 = (TH2D*)fs2.Get(histosyst2[i].c_str());

    for (unsigned int j = 1; j <= histo->GetNbinsX(); j++) {
      for (unsigned int k = 1; k <= histo->GetNbinsY(); k++) {
	if (histo->GetBinContent(j,k) > 0.0001) {

	  cout << std::fixed;
	  int tempdig, precis;
	  if (k == 1) {
	    cout << "\\hline" << endl << endl << endl;
	   //  cout << "$" << std::setprecision(1) << histo->GetXaxis()->GetBinLowEdge(j) 
	//	 << "-" << std::setprecision(1) << histo->GetXaxis()->GetBinUpEdge(j)
	//	 << "$ ";
	  }
	  precis = 1;
	  if (j == 4) precis = 2;
	  // cout << "& $"  << std::setprecision(precis) << histo->GetYaxis()->GetBinLowEdge(k) 
	  //     << "-" << std::setprecision(precis) << histo->GetYaxis()->GetBinUpEdge(k);
	  // Compute Taylor's errors
	  float cont = histo->GetBinContent(j,k);
	  float err = histo->GetBinError(j,k);
	  // float syst = histo2->GetBinContent(j,k)*cont/100.;
	  float temperr, roundedCont, roundedErr, roundedSyst;
	  
	  for (int idig = -5; idig < 6; idig++) {
	    if (err/pow(10,idig) < 10.0) {
	      temperr = err/pow(10,idig);
	      tempdig = idig;
	      break;
	    }
	  }
	  precis = 0;
	  if (int(temperr) == 1 || int(temperr) == 2) precis = 1;
	  int theRounding = precis - tempdig;
	  roundedErr = ceil(err*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  roundedCont = ceil(cont*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  // roundedSyst = ceil(syst*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  /* if (theRounding < 0) theRounding = 0;
	  cout << "$ & $" << std::setprecision(theRounding) << roundedCont 
	       << " \\pm " << std::setprecision(theRounding) << roundedErr;
	  //     << " \\pm " << std::setprecision(theRounding) << roundedSyst;

	  cont = histo3->GetBinContent(j,k);
	  err = histo3->GetBinError(j,k);
	  // syst = histo4->GetBinContent(j,k)*cont/100.;

	  for (int idig = -5; idig < 6; idig++) {
	    if (err/pow(10,idig) < 10.0) {
	      temperr = err/pow(10,idig);
	      tempdig = idig;
	      break;
	    }
	  }
	  precis = 0;
	  if (int(temperr) == 1 || int(temperr) == 2) precis = 1;
	  theRounding = precis - tempdig;
	  roundedErr = ceil(err*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  roundedCont = ceil(cont*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  // roundedSyst = ceil(syst*pow(10.,theRounding) - 0.5)/pow(10.,theRounding);
	  */
	  if (theRounding < 0) theRounding = 0;
	  cout << "$ & $" << std::setprecision(theRounding) << roundedCont 
	       << " \\pm " << std::setprecision(theRounding) << roundedErr
	       // << " \\pm " << std::setprecision(theRounding) << roundedSyst
	       << "$ \\\\" << endl;
	}
      }
    }
  }
  
  return; 
}
