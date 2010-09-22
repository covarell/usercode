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

void computeDiffXsec() {

  static const unsigned int smallBins = 33;
  static const unsigned int largeBins = 15;
  int howManyBins[largeBins] = {2,2,2,2,2,2,1,4,3,3,3,2,2,2,1};
  double binWidths[smallBins] = {1.,1.,1.5,2.,2.,18.,1.5,1.,1.,1.,1.5,2.,20.,
				 0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
				 0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,1.,
				 1.5,2.,20.};
  double binWidths2[largeBins] = {2.,3.5,20.,2.5,2.,3.5,20.,1.25,0.75,0.75,0.75,1.,2.,3.5,20.};
   
  // YIELDS AND ERRORS ("MIXED" POLARIZATION)
  double inclXsec[smallBins] = {35.33,23.29,18.31,7.84,2.86,0.278,
				55.04,37.95,23.62,13.22,6.12,2.21,0.145,
				62.64,146.08,170.64,219.17,203.07,230.44,
				194.60,200.09,157.40,155.64,123.94,118.03,
				107.12,80.03,62.94,38.27,20.22,9.92,3.85,0.239};
  double errstXsec[smallBins] = {2.31,1.38,0.71,0.28,0.13,0.011,
				 4.92,2.33,1.19,0.65,0.24,0.11,0.008,
				 4.37,10.47,10.28,12.18,9.87,11.05,
				 9.69,10.14,31.46,9.66,5.98,5.52,
				 4.80,17.81,2.32,1.19,0.72,0.39,0.17,0.012};
  double errsyXsec[smallBins] = {9.33,5.72,2.93,1.20,0.57,0.044,
				 14.55,9.14,3.56,1.64,0.79,0.30,0.022,
				 15.03,49.32,62.07,73.25,69.14,73.25,
				 50.51,43.32,30.94,28.03,21.03,19.93,
				 16.33,12.23,8.40,5.43,2.81,1.28,0.49,0.031};

  // POLARIZATION FROM FILE
  double nopol[smallBins];
  double trasvHX[smallBins];
  double longHX[smallBins];
  double trasvCS[smallBins];
  double longCS[smallBins];
  double B409[smallBins];
  double B333[smallBins];

  TFile* accfile[3];

  accfile[0] = TFile::Open("accept/0._1.2_acceptance_polarized.root");
  accfile[1] = TFile::Open("accept/1.2_1.6_acceptance_polarized.root");
  accfile[2] = TFile::Open("accept/1.6_2.4_acceptance_polarized.root");

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
  double errsyBfrac[largeBins] = {0.019,0.023,0.014,0.034,0.036,0.054,0.028,0.050,0.024,0.022,0.009,0.012,0.009,0.022,0.019};

  int tsb = 0;
  for (int i=0; i < largeBins; i++) {
    
    double totinclXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclstXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclsyXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double promptXsec[5] = {0.,0.,0.,0.,0.};
    double promptstXsec[5] = {0.,0.,0.,0.,0.};
    double promptsyXsec[5] = {0.,0.,0.,0.,0.};
    double bXsec[2] = {0.,0.};
    double bstXsec[2] = {0.,0.};
    double bsyXsec[2] = {0.,0.};

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
      cout << "Prompt polarization type " << polpr << " : $" << promptXsec[polpr] << " \\pm " << promptstXsec[polpr] << " \\pm " << promptsyXsec[polpr] << "$" << endl;
    }

    for (int polnpr=0; polnpr < 2; polnpr++) {
      bXsec[polnpr] = totinclXsec[polnpr+5]*Bfrac[i];
      bstXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclstXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errstBfrac[i]/Bfrac[i],2));
      bsyXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclsyXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errsyBfrac[i]/Bfrac[i],2));
      cout << "Non-prompt polarization type " << polnpr << " : $" << bXsec[polnpr] << " \\pm " << bstXsec[polnpr] << " \\pm " << bsyXsec[polnpr] << "$" << endl;
    }
      
  }
      
  return;
}
