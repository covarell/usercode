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

  static const unsigned int smallBins = 13;
  static const unsigned int largeBins = 8;
  int howManyBins[largeBins] = {1,2,1,3,2,1,2,1};
  double binWidths[smallBins] = {2.,2.,2.,20.,1.,0.5,0.5,1.,1.,2.,2.,2.,20.};
  double binWidths2[largeBins] = {2.,4.,20.,2.,2.,2.,4.,20.};
   
  // YIELDS AND ERRORS ("MIXED" POLARIZATION)
  // double inclXsec[smallBins] = {34.622,15.951,8.514,0.662,191.050,375.838,344.348,199.419,111.534,55.972,14.289,7.019,0.322};
  // double errstXsec[smallBins] = {2.971,0.936,0.498,0.035,11.849,34.427,38.920,8.197,12.335,2.455,0.690,0.440,0.024};
  // double errsyXsec[smallBins] = {7.575,2.384,1.379,0.091,58.991,95.725,72.548,30.253,18.289,13.247,2.395,2.085,0.050};
  // YIELDS AND ERRORS (NO POLARIZATION)
  double inclXsec[smallBins] = {34.89,16.18,8.49,0.653,214.0,357.0,335.0,197.0,109.0,54.6,14.91,5.88,0.307};
  double errstXsec[smallBins] = {2.50,0.84,0.45,0.031,14.0,34.0,20.0,8.0,4.0,3.1,0.64,0.34,0.024};
  double errsyXsec[smallBins] = {6.00,2.33,1.35,0.097,66.0,90.0,70.0,29.0,17.0,10.0,2.61,1.00,0.048};

  // POLARIZATION FROM FILE
  double nopol[smallBins];
  double trasvHX[smallBins];
  double longHX[smallBins];
  double trasvCS[smallBins];
  double longCS[smallBins];
  double B409[smallBins];
  double B333[smallBins];

  TFile f1("large_bins/acceptance_polarized.root");

  TH2F *Jpsi_non = (TH2F*)f1.Get("Jpsi_non");
  TH2F *Jpsi_t_HXtheta = (TH2F*)f1.Get("Jpsi_t_HXtheta");
  TH2F *Jpsi_l_HXtheta = (TH2F*)f1.Get("Jpsi_l_HXtheta");
  TH2F *Jpsi_t_CStheta = (TH2F*)f1.Get("Jpsi_t_CStheta");
  TH2F *Jpsi_l_CStheta = (TH2F*)f1.Get("Jpsi_l_CStheta");
  TH2F *B_409 = (TH2F*)f1.Get("B_409");
  TH2F *B_333 = (TH2F*)f1.Get("B_333");
  
  int nn = 0;
  for (int ii=1; ii <= Jpsi_non->GetNbinsX(); ii++) {
    for (int jj=1; jj <= Jpsi_non->GetNbinsY(); jj++) {
      if ( ii == 2 || jj > 5 ) {  // barrel region excluded
	nopol[nn]   = Jpsi_non->GetBinContent(ii,jj);
	trasvHX[nn] = Jpsi_t_HXtheta->GetBinContent(ii,jj);	
	longHX[nn]  = Jpsi_l_HXtheta->GetBinContent(ii,jj);
	trasvCS[nn] = Jpsi_t_CStheta->GetBinContent(ii,jj);
	longCS[nn]  = Jpsi_l_CStheta->GetBinContent(ii,jj);
	B409[nn]   = B_409->GetBinContent(ii,jj);
	B333[nn]   = B_333->GetBinContent(ii,jj);
	cout << "nn " << ii << " " << jj << " " << nn << " " << nopol[nn] << " " << trasvHX[nn] << " " << longHX[nn] << " " << trasvCS[nn] << " " << longCS[nn] << " " << B409[nn] << " " << B333[nn] << endl;
	nn++;
      }
    }
  }
  cout << "OK qui " << endl;  

  // B-FRACTION
  double Bfrac[largeBins] = {0.162,0.257,0.369,0.098,0.112,0.165,0.203,0.331};
  double errstBfrac[largeBins] = {0.038,0.022,0.027,0.022,0.013,0.019,0.019,0.039};
  double errsyBfrac[largeBins] = {0.033,0.011,0.014,0.036,0.011,0.010,0.010,0.018};

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
