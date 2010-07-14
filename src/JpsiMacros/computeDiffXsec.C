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
  double binWidths[smallBins] = {2.,2.,4.,20.,1.,0.5,0.5,1.,1.,2.,2.,4.,20.};
  double binWidths2[largeBins] = {2.,6.,20.,2.,2.,2.,6.,20.};
   
  // YIELDS AND ERRORS ("MIXED" POLARIZATION)
  // double inclXsec[smallBins] = {34.622,15.951,8.514,0.662,191.050,375.838,344.348,199.419,111.534,55.972,14.289,7.019,0.322};
  // double errstXsec[smallBins] = {2.971,0.936,0.498,0.035,11.849,34.427,38.920,8.197,12.335,2.455,0.690,0.440,0.024};
  // double errsyXsec[smallBins] = {7.575,2.384,1.379,0.091,58.991,95.725,72.548,30.253,18.289,13.247,2.395,2.085,0.050};
  // YIELDS AND ERRORS (NO POLARIZATION)
  double inclXsec[smallBins] = {34.622,15.951,8.514,0.662,191.050,375.838,344.348,199.419,111.534,55.972,14.289,7.019,0.322};
  double errstXsec[smallBins] = {2.971,0.936,0.498,0.035,11.849,34.427,38.920,8.197,12.335,2.455,0.690,0.440,0.024};
  double errsyXsec[smallBins] = {7.575,2.384,1.379,0.091,58.991,95.725,72.548,30.253,18.289,13.247,2.395,2.085,0.050};

  // POLARIZATION FROM FILE
  double nopol[smallBins];
  double trasvHX[smallBins];
  double longHX[smallBins];
  double trasvCS[smallBins];
  double longCS[smallBins];
  double B409[smallBins];

  TFile f1("large_bins/acceptance_polarized.root");

  TH2F *Jpsi_non = (TH2F*)f1.Get("Jpsi_non");
  TH2F *Jpsi_t_HXtheta = (TH2F*)f1.Get("Jpsi_t_HXtheta");
  TH2F *Jpsi_l_HXtheta = (TH2F*)f1.Get("Jpsi_l_HXtheta");
  TH2F *Jpsi_t_CStheta = (TH2F*)f1.Get("Jpsi_t_CStheta");
  TH2F *Jpsi_l_CStheta = (TH2F*)f1.Get("Jpsi_l_CStheta");
  TH2F *B_409 = (TH2F*)f1.Get("B_409");
  
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
	cout << "nn " << ii << " " << jj << " " << nn << " " << nopol[nn] << " " << trasvHX[nn] << " " << longHX[nn] << " " << trasvCS[nn] << " " << longCS[nn] << " " << B409[nn] << endl;
	nn++;
      }
    }
  }
  cout << "OK qui " << endl;  

  // B-FRACTION
  double Bfrac[largeBins] = {0.158,0.255,0.375,0.085,0.115,0.189,0.207,0.346};
  double errstBfrac[largeBins] = {0.042,0.024,0.029,0.023,0.015,0.021,0.021,0.043};
  // to be updated
  double errsyBfrac[largeBins] = {0.018,0.009,0.012,0.027,0.016,0.011,0.017,0.010};

  int tsb = 0;
  for (int i=0; i < largeBins; i++) {
    
    double totinclXsec[6] = {0.,0.,0.,0.,0.,0.};
    double totinclstXsec[6] = {0.,0.,0.,0.,0.,0.};
    double totinclsyXsec[6] = {0.,0.,0.,0.,0.,0.};
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
      tsb++;
    }

    cout << endl;
    cout << "Bin number " << i << endl << endl;    

    for (int pol=0; pol < 6; pol++) {
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

    bXsec[0] = totinclXsec[0]*Bfrac[i];
    bstXsec[0] = bXsec[0]*sqrt(pow(totinclstXsec[0]/totinclXsec[0],2) + pow(errstBfrac[i]/Bfrac[i],2));
    bsyXsec[0] = bXsec[0]*sqrt(pow(totinclsyXsec[0]/totinclXsec[0],2) + pow(errsyBfrac[i]/Bfrac[i],2));
    bXsec[1] = totinclXsec[5]*Bfrac[i];
    bstXsec[1] = bXsec[1]*sqrt(pow(totinclstXsec[5]/totinclXsec[5],2) + pow(errstBfrac[i]/Bfrac[i],2));
    bsyXsec[1] = bXsec[1]*sqrt(pow(totinclsyXsec[5]/totinclXsec[5],2) + pow(errsyBfrac[i]/Bfrac[i],2));

    for (int polnpr=0; polnpr < 2; polnpr++) {
      cout << "Non-prompt polarization type " << polnpr << " : $" << bXsec[polnpr] << " \\pm " << bstXsec[polnpr] << " \\pm " << bsyXsec[polnpr] << "$" << endl;
    }
      
  }
      
  return;
}
