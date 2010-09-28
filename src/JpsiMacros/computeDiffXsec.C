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

void computeDiffXsec(bool isDoubleDiff = false) { 

  static const unsigned int smallBins = 33;
  static const unsigned int largeBins = 15;
  int howManyBins[largeBins] = {2,2,2,2,2,2,1,4,3,3,3,2,2,2,1};
  double binWidths[smallBins] = {1.,1.,1.5,2.,2.,18.,1.5,1.,1.,1.,1.5,2.,20.,
				 0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
				 0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,1.,
				 1.5,2.,20.};
  double binWidths2[largeBins] = {2.,3.5,20.,2.5,2.,3.5,20.,1.25,0.75,0.75,0.75,1.,2.,3.5,20.};

  double ybinWidths[3] = {2.4,0.8,1.6};
  if (!isDoubleDiff) {
    for (int k=0; k < 3; k++) { ybinWidths[k] = 1.0; }
  }
   
  // YIELDS AND ERRORS ("MIXED" POLARIZATION)
  double inclXsec[smallBins] = {35.2,23.3,18.3,7.75,2.83,0.278,
				55,36.9,22.8,13.2,6.11,2.21,0.145,
				58.9,133,164,195,204,212,
				195,200,154,154,125,118,
				107,79.4,63.5,39.2,20.1,9.91,3.85,0.239};
  double errstXsec[smallBins] = {2.3,1.4,0.71,0.27,0.13,0.011,
				 5,2.1,1.1,0.65,0.24,0.11,0.008,
				 3.5,7.2,7.9,8.5,8.9,8.5,
				 9.9,9.8,6.8,12,5.9,5.6,
				 5.1,2.8,2.3,1.2,0.7,0.39,0.17,0.012};
  double errsyXsec[smallBins] = {5.5,3.3,2.4,0.9,0.33,0.031,
				 9.4,4.9,3.0,1.6,0.7,0.25,0.017,
				 8.5,22,24,31,31,32,
				 25,27,21,19,16,14,
				 13,10,7.4,5,2.6,1.2,0.45,0.031};

  // POLARIZATION FROM FILE
  double nopol[smallBins];
  double trasvHX[smallBins];
  double longHX[smallBins];
  double trasvCS[smallBins];
  double longCS[smallBins];
  double B409[smallBins];
  double B333[smallBins];

  TFile* accfile[3];

  accfile[0] = TFile::Open("accept/0.0_1.2_Acc_acceptance_polarized.root");
  accfile[1] = TFile::Open("accept/1.2_1.6_Acc_acceptance_polarized.root");
  accfile[2] = TFile::Open("accept/1.6_2.4_Acc_acceptance_polarized.root");

  TH2F *Jpsi_non[3];
  TH2F *Jpsi_t_HXtheta[3];
  TH2F *Jpsi_l_HXtheta[3];
  TH2F *Jpsi_t_CStheta[3];
  TH2F *Jpsi_l_CStheta[3];
  TH2F *B_409[3];
  TH2F *B_333[3];

  for (int i=0; i < 3; i++) {
    Jpsi_non[i] = (TH2F*)(accfile[i]->Get("Jpsi_non"));
    Jpsi_t_HXtheta[i] = (TH2F*)(accfile[i]->Get("Jpsi_t_HXtheta"));
    Jpsi_l_HXtheta[i] = (TH2F*)(accfile[i]->Get("Jpsi_l_HXtheta"));
    Jpsi_t_CStheta[i] = (TH2F*)(accfile[i]->Get("Jpsi_t_CStheta"));
    Jpsi_l_CStheta[i] = (TH2F*)(accfile[i]->Get("Jpsi_l_CStheta"));
    B_409[i] = (TH2F*)(accfile[i]->Get("B_409"));
    B_333[i] = (TH2F*)(accfile[i]->Get("B_333"));
  }
  
  int nn = 0;
  for (int ii=1; ii <= 3; ii++) {
    for (int jj=1; jj <= Jpsi_non[ii-1]->GetNbinsY(); jj++) {
      nopol[nn]   = Jpsi_non[ii-1]->GetBinContent(ii,jj);
      trasvHX[nn] = Jpsi_t_HXtheta[ii-1]->GetBinContent(ii,jj);	
      longHX[nn]  = Jpsi_l_HXtheta[ii-1]->GetBinContent(ii,jj);
      trasvCS[nn] = Jpsi_t_CStheta[ii-1]->GetBinContent(ii,jj);
      longCS[nn]  = Jpsi_l_CStheta[ii-1]->GetBinContent(ii,jj);
      B409[nn]   = B_409[ii-1]->GetBinContent(ii,jj);
      B333[nn]   = B_333[ii-1]->GetBinContent(ii,jj);
      cout << "nn " << ii << " " << jj << " " << nn << " " << nopol[nn] << " " << trasvHX[nn] << " " << longHX[nn] << " " << trasvCS[nn] << " " << longCS[nn] << " " << B409[nn] << " " << B333[nn] << endl;
      nn++;
    }
  }
  cout << "OK qui " << endl;  
  
  // B-FRACTION
  double Bfrac[largeBins] = {0.181,0.261,0.405,0.147,0.178,0.204,0.363,0.057,0.087,0.113,0.138,0.160,0.173,0.233,0.377};
  double errstBfrac[largeBins] = {0.024,0.015,0.018,0.021,0.017,0.017,0.031,0.021,0.014,0.013,0.014,0.014,0.012,0.016,0.031};
  double errsyBfrac[largeBins] = {0.014,0.014,0.005,0.027,0.019,0.003,0.015,
0.042,0.022,0.020,0.009,0.012,0.010,0.012,0.005};

  int tsb = 0;

  double totpr, errsttotpr, errsytotpr;
  double totb, errsttotb, errsytotb;

  for (int i=0; i < largeBins; i++) {
    
    double totinclXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclstXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclsyXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double promptXsec[5] = {0.,0.,0.,0.,0.};
    double promptstXsec[5] = {0.,0.,0.,0.,0.};
    double promptsyXsec[5] = {0.,0.,0.,0.,0.};
    double prompttoXsec[5] = {0.,0.,0.,0.,0.};
    double bXsec[2] = {0.,0.};
    double bstXsec[2] = {0.,0.};
    double bsyXsec[2] = {0.,0.};
    double btoXsec[2] = {0.,0.};

    for (int k=0; k < howManyBins[i]; k++) {
      totinclXsec[0] += binWidths[tsb]*inclXsec[tsb];
      totinclstXsec[0] += pow(binWidths[tsb]*errstXsec[tsb],2);
      totinclsyXsec[0] += pow(binWidths[tsb]*errsyXsec[tsb],2);
      totinclXsec[4] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/trasvHX[tsb];
      totinclstXsec[4] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/trasvHX[tsb],2);
      totinclsyXsec[4] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/trasvHX[tsb],2);
      totinclXsec[2] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/trasvCS[tsb];
      totinclstXsec[2] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/trasvCS[tsb],2);
      totinclsyXsec[2] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/trasvCS[tsb],2);
      totinclXsec[3] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/longHX[tsb];
      totinclstXsec[3] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/longHX[tsb],2);
      totinclsyXsec[3] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/longHX[tsb],2);
      totinclXsec[1] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/longCS[tsb];
      totinclstXsec[1] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/longCS[tsb],2);
      totinclsyXsec[1] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/longCS[tsb],2);
      totinclXsec[5] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/B409[tsb];
      totinclstXsec[5] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/B409[tsb],2);
      totinclsyXsec[5] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/B409[tsb],2);
      totinclXsec[6] += binWidths[tsb]*inclXsec[tsb]*nopol[tsb]/B333[tsb];
      totinclstXsec[6] += pow(binWidths[tsb]*errstXsec[tsb]*nopol[tsb]/B333[tsb],2);
      totinclsyXsec[6] += pow(binWidths[tsb]*errsyXsec[tsb]*nopol[tsb]/B333[tsb],2);
      tsb++;
    }

    cout << endl;
    cout << "Bin number " << i << endl << endl;    

    for (int pol=0; pol < 7; pol++) {
      totinclstXsec[pol] = sqrt(totinclstXsec[pol]);
      totinclsyXsec[pol] = sqrt(totinclsyXsec[pol]);
      totinclXsec[pol] /= binWidths2[i];
      totinclstXsec[pol] /= binWidths2[i];
      totinclsyXsec[pol] /= binWidths2[i];
    }

    for (int polpr=0; polpr < 5; polpr++) {
      promptXsec[polpr] = totinclXsec[polpr]*(1-Bfrac[i]);
      promptstXsec[polpr] = promptXsec[polpr]*sqrt(pow(totinclstXsec[polpr]/totinclXsec[polpr],2) + pow(errstBfrac[i]/(1-Bfrac[i]),2));
      promptsyXsec[polpr] = promptXsec[polpr]*sqrt(pow(totinclsyXsec[polpr]/totinclXsec[polpr],2) + pow(errsyBfrac[i]/(1-Bfrac[i]),2));
      prompttoXsec[polpr] = sqrt(pow(promptstXsec[polpr],2) + pow(promptsyXsec[polpr],2));
      // cout << "Prompt polarization type " << polpr << " : $" << promptXsec[polpr] << " \\pm " << promptstXsec[polpr] << " \\pm " << promptsyXsec[polpr] << "$" << endl;
    }

    for (int polnpr=0; polnpr < 2; polnpr++) {
      bXsec[polnpr] = totinclXsec[polnpr+5]*Bfrac[i];
      bstXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclstXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errstBfrac[i]/Bfrac[i],2));
      bsyXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclsyXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errsyBfrac[i]/Bfrac[i],2));
      btoXsec[polnpr] = sqrt(pow(bstXsec[polnpr],2) + pow(bsyXsec[polnpr],2));
      // cout << "Non-prompt polarization type " << polnpr << " : $" << bXsec[polnpr] << " \\pm " << bstXsec[polnpr] << " \\pm " << bsyXsec[polnpr] << "$" << endl;
    }

    cout << "Prompt : $" << promptXsec[0] << " \\pm " << promptstXsec[0] << " \\pm " << promptsyXsec[0] << "$ & $" << promptXsec[1] << " \\pm " << prompttoXsec[1] << "$ & $" << promptXsec[2] << " \\pm " << prompttoXsec[2] << "$ & $" << promptXsec[3] << " \\pm " << prompttoXsec[3] << "$ & $" << promptXsec[4] << " \\pm " << prompttoXsec[4] << "$" << endl;
    cout << "Non-prompt : $" << bXsec[0] << " \\pm " << bstXsec[0] << " \\pm " << bsyXsec[0] << "$ & $" << bXsec[1] << " \\pm " << btoXsec[1] << "$" << endl;
    
    if (i != 3 && (i < 7 || i > 11)) {

      int ybin = 2;
      if (i < 3) ybin = 0;
      else if (i < 7) ybin = 1;

      totpr += promptXsec[0]*binWidths2[i]*ybinWidths[ybin];  
      totb += bXsec[0]*binWidths2[i]*ybinWidths[ybin];
      
      errsttotpr += pow(promptstXsec[0]*binWidths2[i]*ybinWidths[ybin],2);
      errsttotb += pow(bstXsec[0]*binWidths2[i]*ybinWidths[ybin],2);
      
      errsytotpr += pow(promptsyXsec[0]*binWidths2[i]*ybinWidths[ybin],2);
      errsytotb += pow(bsyXsec[0]*binWidths2[i]*ybinWidths[ybin],2);
    } 
  }  

  cout << endl;
  cout << "TOTAL PROMPT: $" << totpr << " \\pm " << sqrt(errsttotpr) << " \\pm " << sqrt(errsytotpr) << "$" << endl;
  cout << "TOTAL NON-PROMPT: $" << totb << " \\pm " << sqrt(errsttotb) << " \\pm " << sqrt(errsytotb) << "$" << endl;

  return;
}
