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

void computeDiffXsec(bool isDoubleDiff = false, bool plotDellaBfrac = false) { 

  static const unsigned int smallBins = 33;
  static const unsigned int largeBins = 15;
  static const unsigned int largeBins0 = 3;
  static const unsigned int largeBins1 = 4;
  static const unsigned int largeBins2 = 8;
  static const unsigned int highPtLargeBins = 3;

  int howManyBins[largeBins] = {2,2,2,2,2,2,1,4,3,3,3,2,2,2,1};
  double binWidths[smallBins] = {1.,1.,1.5,2.,2.,18.,1.5,1.,1.,1.,1.5,2.,20.,
				 0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
				 0.25,0.25,0.25,0.25,0.25,0.5,0.5,1.,1.,
				 1.5,2.,20.};
  double binWidths2[largeBins] = {2.,3.5,20.,2.5,2.,3.5,20.,1.25,0.75,
				  0.75,0.75,1.,2.,3.5,20.};
  double binAverages0[largeBins0] = {5.46,8.14,13.50};
  double binAverages1[largeBins1] = {3.27,5.48,7.89,12.96};
  double binAverages2[largeBins2] = {0.79,1.60,2.35,3.10,3.96,5.35,7.86,13.11};
  double binRMS0[largeBins0] = {0.56,0.97,3.53};
  double binRMS1[largeBins1] = {0.75,0.55,0.93,3.06};
  double binRMS2[largeBins2] = {0.29,0.21,0.22,0.21,0.29,0.57,0.97,3.23};
  double Bfracs0[largeBins0] = {0.,0.,0.};
  double Bfracs1[largeBins1] = {0.,0.,0.,0.};
  double Bfracs2[largeBins2] = {0.,0.,0.,0.,0.,0.,0.,0.};
  double Bfracerrs0[largeBins0] = {0.,0.,0.};
  double Bfracerrs1[largeBins1] = {0.,0.,0.,0.};
  double Bfracerrs2[largeBins2] = {0.,0.,0.,0.,0.,0.,0.,0.};

  double ybinWidths[3] = {2.4,0.8,1.6};
  if (!isDoubleDiff) {
    for (int k=0; k < 3; k++) { ybinWidths[k] = 1.0; }
  }
   
  // YIELDS AND ERRORS ("MIXED" POLARIZATION)
  double inclXsec[smallBins] = {14.68,9.70,7.62,3.23,1.18,0.116,
				68.77,46.08,28.55,16.51,7.64,2.77,0.182,
				36.81,83.21,102.3,121.9,127.7,132.5,
				121.9,125.1,96.33,96.41,77.86,73.72,
				66.73,49.64,39.68,24.50,12.58,6.20,2.41,0.149};
  double errstXsec[smallBins] = {0.95,0.58,0.30,0.11,0.54,0.005,
				 6.30,2.66,1.33,0.81,0.30,0.14,0.010,
				 2.16,4.48,4.96,5.29,5.56,5.32,
				 6.19,6.13,4.22,7.72,3.72,3.48,
				 3.19,1.73,1.43,0.72,0.44,0.24,0.11,0.008};
  double errsyuncXsec[smallBins] = {1.58,0.66,0.39,0.096,0.033,0.002,
				 10.05,3.60,1.49,0.56,0.20,0.073,0.005,
				 3.75,9.46,11.54,13.77,12.55,13.39,
				 9.96,10.89,8.06,6.28,5.47,4.76,
				 3.69,2.75,1.83,1.09,0.61,0.263,0.066,0.007};
  double errsycXsec[smallBins] = {1.73,1.15,0.90,0.37,0.137,0.013,
				 8.37,5.42,3.56,1.88,0.85,0.31,0.020,
				 4.64,12.0,12.4,15.9,17.5,17.3,
				 14.8,15.2,11.5,11.4,9.24,8.83,
				 8.01,6.5,4.7,3.1,1.54,0.70,0.276,0.018};

  // POLARIZATION FROM FILE
  double mixpol[smallBins];
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

  TFile* accMixfile[3];

  accMixfile[0] = TFile::Open("accept/shuang_Acc_Dimuon_mix_centralvalue_0.0.root");
  accMixfile[1] = TFile::Open("accept/shuang_Acc_Dimuon_mix_centralvalue_1.2.root");
  accMixfile[2] = TFile::Open("accept/shuang_Acc_Dimuon_mix_centralvalue_1.6.root");

  TH2F *Jpsi_mix[3];
  TH2F *Jpsi_non[3];
  TH2F *Jpsi_t_HXtheta[3];
  TH2F *Jpsi_l_HXtheta[3];
  TH2F *Jpsi_t_CStheta[3];
  TH2F *Jpsi_l_CStheta[3];
  TH2F *B_409[3];
  TH2F *B_333[3];

  for (int i=0; i < 3; i++) {
    Jpsi_mix[i] = (TH2F*)(accMixfile[i]->Get("Jpsi_non"));
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
      mixpol[nn]   = Jpsi_mix[ii-1]->GetBinContent(ii,jj);
      nopol[nn]   = Jpsi_non[ii-1]->GetBinContent(ii,jj);
      trasvHX[nn] = Jpsi_t_HXtheta[ii-1]->GetBinContent(ii,jj);	
      longHX[nn]  = Jpsi_l_HXtheta[ii-1]->GetBinContent(ii,jj);
      trasvCS[nn] = Jpsi_t_CStheta[ii-1]->GetBinContent(ii,jj);
      longCS[nn]  = Jpsi_l_CStheta[ii-1]->GetBinContent(ii,jj);
      B409[nn]   = B_409[ii-1]->GetBinContent(ii,jj);
      B333[nn]   = B_333[ii-1]->GetBinContent(ii,jj);
      cout << "nn " << ii << " " << jj << " " << nn << " " << mixpol[nn] << " " << nopol[nn] << " " << trasvHX[nn] << " " << longHX[nn] << " " << trasvCS[nn] << " " << longCS[nn] << " " << B409[nn] << " " << B333[nn] << endl;
      nn++;
    }
  }
  cout << "OK qui " << endl;  
  
  // B-FRACTION
  double Bfrac[largeBins] = {0.181,0.261,0.405,0.147,0.178,0.204,0.363,0.057,0.087,0.113,0.138,0.160,0.173,0.233,0.377};
  double errstBfrac[largeBins] = {0.024,0.015,0.018,0.021,0.017,0.017,0.031,0.021,0.014,0.013,0.014,0.014,0.012,0.016,0.031};
  double errsyBfrac[largeBins] = {0.015,0.014,0.005,0.028,0.019,0.004,0.016,
0.042,0.022,0.020,0.009,0.013,0.012,0.012,0.008};

  int tsb = 0;

  double promptXsecint24[highPtLargeBins] = {0.,0.,0.};
  double promptstXsecint24[highPtLargeBins] = {0.,0.,0.};
  double promptsyXsecint24[highPtLargeBins] = {0.,0.,0.};
  double promptsycXsecint24[highPtLargeBins] = {0.,0.,0.};
  double promptsyuncXsecint24[highPtLargeBins] = {0.,0.,0.};
  double prompttoXsecint24[highPtLargeBins] = {0.,0.,0.};
  double bXsecint24[highPtLargeBins] = {0.,0.,0.};
  double bstXsecint24[highPtLargeBins] = {0.,0.,0.};
  double bsyXsecint24[highPtLargeBins] = {0.,0.,0.};
  double bsycXsecint24[highPtLargeBins] = {0.,0.,0.};
  double bsyuncXsecint24[highPtLargeBins] = {0.,0.,0.};
  double btoXsecint24[highPtLargeBins] = {0.,0.,0.};
  double totpr, errsttotpr, errsytotpr; 
  double errsyctotpr, errsyunctotpr;
  double totb, errsttotb, errsytotb; 
  double errsyctotb, errsyunctotb;

  unsigned int thePtBin = 0;

  for (int i=0; i < largeBins; i++) {
    
    double totinclXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclstXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclsycXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double totinclsyuncXsec[7] = {0.,0.,0.,0.,0.,0.,0.};
    double promptXsec[5] = {0.,0.,0.,0.,0.};
    double promptstXsec[5] = {0.,0.,0.,0.,0.};
    double promptsyXsec[5] = {0.,0.,0.,0.,0.};
    double promptsycXsec[5] = {0.,0.,0.,0.,0.};
    double promptsyuncXsec[5] = {0.,0.,0.,0.,0.};
    double prompttoXsec[5] = {0.,0.,0.,0.,0.};
    double bXsec[2] = {0.,0.};
    double bstXsec[2] = {0.,0.};
    double bsyXsec[2] = {0.,0.};
    double bsycXsec[2] = {0.,0.};
    double bsyuncXsec[2] = {0.,0.};
    double btoXsec[2] = {0.,0.};

    for (int k=0; k < howManyBins[i]; k++) {
      totinclXsec[0] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/nopol[tsb];
      totinclstXsec[0] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/nopol[tsb],2);
      totinclsycXsec[0] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/nopol[tsb];;
      totinclsyuncXsec[0] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/nopol[tsb],2);

      totinclXsec[4] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/trasvHX[tsb];
      totinclstXsec[4] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/trasvHX[tsb],2);
      totinclsycXsec[4] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/trasvHX[tsb];
      totinclsyuncXsec[4] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/trasvHX[tsb],2);

      totinclXsec[2] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/trasvCS[tsb];
      totinclstXsec[2] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/trasvCS[tsb],2);
      totinclsycXsec[2] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/trasvCS[tsb];
      totinclsyuncXsec[2] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/trasvCS[tsb],2);

      totinclXsec[3] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/longHX[tsb];
      totinclstXsec[3] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/longHX[tsb],2);
      totinclsycXsec[3] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/longHX[tsb];
      totinclsyuncXsec[3] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/longHX[tsb],2);

      totinclXsec[1] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/longCS[tsb];
      totinclstXsec[1] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/longCS[tsb],2);
      totinclsycXsec[1] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/longCS[tsb];
      totinclsyuncXsec[1] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/longCS[tsb],2);

      totinclXsec[5] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/B409[tsb];
      totinclstXsec[5] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/B409[tsb],2);
      totinclsycXsec[5] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/B409[tsb];
      totinclsyuncXsec[5] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/B409[tsb],2); 

      totinclXsec[6] += binWidths[tsb]*inclXsec[tsb]*mixpol[tsb]/B333[tsb];
      totinclstXsec[6] += pow(binWidths[tsb]*errstXsec[tsb]*mixpol[tsb]/B333[tsb],2);
      totinclsycXsec[6] += binWidths[tsb]*errsycXsec[tsb]*mixpol[tsb]/B333[tsb];
      totinclsyuncXsec[6] += pow(binWidths[tsb]*errsyuncXsec[tsb]*mixpol[tsb]/B333[tsb],2);
      
      tsb++;
    }

    cout << endl;
    cout << "Bin number " << i << endl << endl;    

    for (int pol=0; pol < 7; pol++) {
      totinclstXsec[pol] = sqrt(totinclstXsec[pol]);
      totinclsyuncXsec[pol] = sqrt(totinclsyuncXsec[pol]);
      totinclXsec[pol] /= binWidths2[i];
      totinclstXsec[pol] /= binWidths2[i];
      totinclsycXsec[pol] /= binWidths2[i];
      totinclsyuncXsec[pol] /= binWidths2[i];
    }

    for (int polpr=0; polpr < 5; polpr++) {
      promptXsec[polpr] = totinclXsec[polpr]*(1-Bfrac[i]);
      promptstXsec[polpr] = promptXsec[polpr]*sqrt(pow(totinclstXsec[polpr]/totinclXsec[polpr],2) + pow(errstBfrac[i]/(1-Bfrac[i]),2));
      promptsycXsec[polpr] = promptXsec[polpr]*sqrt(pow(totinclsycXsec[polpr]/totinclXsec[polpr],2) + pow(errsyBfrac[i]/(1-Bfrac[i]),2));
      promptsyuncXsec[polpr] = promptXsec[polpr]*sqrt(pow(totinclsyuncXsec[polpr]/totinclXsec[polpr],2) + pow(errsyBfrac[i]/(1-Bfrac[i]),2));
      promptsyXsec[polpr] = sqrt(pow(promptsycXsec[polpr],2) + pow(promptsyuncXsec[polpr],2));
      prompttoXsec[polpr] = sqrt(pow(promptstXsec[polpr],2) + pow(promptsyXsec[polpr],2));
      // cout << "Prompt polarization type " << polpr << " : $" << promptXsec[polpr] << " \\pm " << promptstXsec[polpr] << " \\pm " << promptsyXsec[polpr] << "$" << endl;
    }

    for (int polnpr=0; polnpr < 2; polnpr++) {
      bXsec[polnpr] = totinclXsec[polnpr+5]*Bfrac[i];
      bstXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclstXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errstBfrac[i]/Bfrac[i],2));
      bsycXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclsycXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errsyBfrac[i]/Bfrac[i],2));
      bsyuncXsec[polnpr] = bXsec[polnpr]*sqrt(pow(totinclsyuncXsec[polnpr+5]/totinclXsec[polnpr+5],2) + pow(errsyBfrac[i]/Bfrac[i],2));
      bsyXsec[polnpr] = sqrt(pow(bsycXsec[polnpr],2) + pow(bsyuncXsec[polnpr],2));
      btoXsec[polnpr] = sqrt(pow(bstXsec[polnpr],2) + pow(bsyXsec[polnpr],2));
      // cout << "Non-prompt polarization type " << polnpr << " : $" << bXsec[polnpr] << " \\pm " << bstXsec[polnpr] << " \\pm " << bsyXsec[polnpr] << "$" << endl;
    }

    cout << "Prompt : $" << promptXsec[0] << " \\pm " << prompttoXsec[0] << " (\\pm " << promptstXsec[0] << "_{\\mathrm{stat.}} \\pm " << promptsyXsec[0] << "_{\\mathrm{syst.}})$ & $" << promptXsec[1] << " \\pm " << prompttoXsec[1] << "$ & $" << promptXsec[2] << " \\pm " << prompttoXsec[2] << "$ & $" << promptXsec[3] << " \\pm " << prompttoXsec[3] << "$ & $" << promptXsec[4] << " \\pm " << prompttoXsec[4] << "$" << endl;
    cout << "Non-prompt : $" << bXsec[0] << " \\pm " << btoXsec[0] << " (\\pm " << bstXsec[0] << "_{\\mathrm{stat.}} \\pm " << bsyXsec[0] << "_{\\mathrm{syst.}})$" << endl;
    double BfracCorr = bXsec[0]/(promptXsec[0]+bXsec[0]);
    double totErr = sqrt(pow(errstBfrac[i]*(BfracCorr/Bfrac[i]),2) + pow(errsyBfrac[i]*(BfracCorr/Bfrac[i]),2));
    if (plotDellaBfrac) {
      cout << "Corrected B-fraction : $" << BfracCorr << " \\pm " << errstBfrac[i]*(BfracCorr/Bfrac[i]) << "_{\\mathrm{stat.}} \\pm " << errsyBfrac[i]*(BfracCorr/Bfrac[i]) << "_{\\mathrm{syst.}})$" << endl; 
    }
    
    if (i < 3) {

      thePtBin = i;
      promptXsecint24[thePtBin] += promptXsec[0]*ybinWidths[0];  
      bXsecint24[thePtBin] += bXsec[0]*ybinWidths[0];
      promptstXsecint24[thePtBin] += pow(promptstXsec[0]*ybinWidths[0],2);
      bstXsecint24[thePtBin]+= pow(bstXsec[0]*ybinWidths[0],2);
      promptsyuncXsecint24[thePtBin] += pow(promptsyuncXsec[0]*ybinWidths[0],2);
      bsyuncXsecint24[thePtBin] += pow(bsyuncXsec[0]*ybinWidths[0],2);
      promptsycXsecint24[thePtBin] += promptsycXsec[0]*ybinWidths[0];
      bsycXsecint24[thePtBin] += bsycXsec[0]*ybinWidths[0];

      totpr += promptXsec[0]*binWidths2[i]*ybinWidths[0];  
      totb += bXsec[0]*binWidths2[i]*ybinWidths[0];
      errsttotpr += pow(promptstXsec[0]*binWidths2[i]*ybinWidths[0],2);
      errsttotb += pow(bstXsec[0]*binWidths2[i]*ybinWidths[0],2);
      errsyunctotpr += pow(promptsyuncXsec[0]*binWidths2[i]*ybinWidths[0],2);
      errsyunctotb += pow(bsyuncXsec[0]*binWidths2[i]*ybinWidths[0],2);
      errsyctotpr += promptsycXsec[0]*binWidths2[i]*ybinWidths[0];
      errsyctotb += bsycXsec[0]*binWidths2[i]*ybinWidths[0];

      if (plotDellaBfrac) {
	Bfracs0[i] = BfracCorr;
        Bfracerrs0[i] = totErr;
      }

    } 

    if (i > 2 && i < 7 && plotDellaBfrac) {
      Bfracs1[i-3] = BfracCorr;
      Bfracerrs1[i-3] = totErr;
    }

    if (i > 3 && i < 7) {

      thePtBin = i-4;
      promptXsecint24[thePtBin] += promptXsec[0]*ybinWidths[1];  
      bXsecint24[thePtBin] += bXsec[0]*ybinWidths[1];
      promptstXsecint24[thePtBin] += pow(promptstXsec[0]*ybinWidths[1],2);
      bstXsecint24[thePtBin]+= pow(bstXsec[0]*ybinWidths[1],2);
      promptsyuncXsecint24[thePtBin] += pow(promptsyuncXsec[0]*ybinWidths[1],2);
      bsyuncXsecint24[thePtBin] += pow(bsyuncXsec[0]*ybinWidths[1],2);
      promptsycXsecint24[thePtBin] += promptsycXsec[0]*ybinWidths[1];
      bsycXsecint24[thePtBin] += bsycXsec[0]*ybinWidths[1];

      totpr += promptXsec[0]*binWidths2[i]*ybinWidths[1];  
      totb += bXsec[0]*binWidths2[i]*ybinWidths[1];
      errsttotpr += pow(promptstXsec[0]*binWidths2[i]*ybinWidths[1],2);
      errsttotb += pow(bstXsec[0]*binWidths2[i]*ybinWidths[1],2);
      errsyunctotpr += pow(promptsyuncXsec[0]*binWidths2[i]*ybinWidths[1],2);
      errsyunctotb += pow(bsyuncXsec[0]*binWidths2[i]*ybinWidths[1],2);
      errsyctotpr += promptsycXsec[0]*binWidths2[i]*ybinWidths[1];
      errsyctotb += bsycXsec[0]*binWidths2[i]*ybinWidths[1];

    } 

    if (i > 6 && plotDellaBfrac) {
      Bfracs2[i-7] = BfracCorr;
      Bfracerrs2[i-7] = totErr;
    }

    if (i > 11) {

      thePtBin = i-12;
      promptXsecint24[thePtBin] += promptXsec[0]*ybinWidths[2];  
      bXsecint24[thePtBin] += bXsec[0]*ybinWidths[2];
      promptstXsecint24[thePtBin] += pow(promptstXsec[0]*ybinWidths[2],2);
      bstXsecint24[thePtBin]+= pow(bstXsec[0]*ybinWidths[2],2);
      promptsyuncXsecint24[thePtBin] += pow(promptsyuncXsec[0]*ybinWidths[2],2);
      bsyuncXsecint24[thePtBin] += pow(bsyuncXsec[0]*ybinWidths[2],2);
      promptsycXsecint24[thePtBin] += promptsycXsec[0]*ybinWidths[2];
      bsycXsecint24[thePtBin] += bsycXsec[0]*ybinWidths[2];

      totpr += promptXsec[0]*binWidths2[i]*ybinWidths[2];  
      totb += bXsec[0]*binWidths2[i]*ybinWidths[2];
      errsttotpr += pow(promptstXsec[0]*binWidths2[i]*ybinWidths[2],2);
      errsttotb += pow(bstXsec[0]*binWidths2[i]*ybinWidths[2],2);
      errsyunctotpr += pow(promptsyuncXsec[0]*binWidths2[i]*ybinWidths[2],2);
      errsyunctotb += pow(bsyuncXsec[0]*binWidths2[i]*ybinWidths[2],2);
      errsyctotpr += promptsycXsec[0]*binWidths2[i]*ybinWidths[2];
      errsyctotb += bsycXsec[0]*binWidths2[i]*ybinWidths[2];

    } 
  }  

  for (int w=0; w < highPtLargeBins; w++) {

    cout << endl;
    cout << "Combined bin number " << w << endl << endl;    

    promptstXsecint24[w] = sqrt(promptstXsecint24[w]);
    promptsyuncXsecint24[w] = sqrt(promptsyuncXsecint24[w]);
    promptXsecint24[w] /= 4.8;
    promptstXsecint24[w] /= 4.8;
    promptsycXsecint24[w] /= 4.8;
    promptsyuncXsecint24[w] /= 4.8;

    bstXsecint24[w] = sqrt(bstXsecint24[w]);
    bsyuncXsecint24[w] = sqrt(bsyuncXsecint24[w]);
    bXsecint24[w] /= 4.8;
    bstXsecint24[w] /= 4.8;
    bsycXsecint24[w] /= 4.8;
    bsyuncXsecint24[w] /= 4.8;

    promptsyXsecint24[w] = sqrt(pow(promptsycXsecint24[w],2) + pow(promptsyuncXsecint24[w],2));
    bsyXsecint24[w] = sqrt(pow(bsycXsecint24[w],2) + pow(bsyuncXsecint24[w],2));
    prompttoXsecint24[w] = sqrt(pow(promptstXsecint24[w],2) + pow(promptsyXsecint24[w],2));
    btoXsecint24[w] = sqrt(pow(bstXsecint24[w],2) + pow(bsyXsecint24[w],2));

    cout << "Prompt : $" << promptXsecint24[w] << " \\pm " << prompttoXsecint24[w] << " (\\pm " << promptstXsecint24[w] << "_{\\mathrm{stat.}} \\pm " << promptsyXsecint24[w] << "_{\\mathrm{syst.}})$" << endl;
    cout << "Non-prompt : $" << bXsecint24[w] << " \\pm " << btoXsecint24[w] << " (\\pm " << bstXsecint24[w] << "_{\\mathrm{stat.}} \\pm " << bsyXsecint24[w] << "_{\\mathrm{syst.}})$" << endl;
  }

  errsytotpr = errsyunctotpr + pow(errsyctotpr,2);
  errsytotb = errsyunctotb + pow(errsyctotb,2);
  //  cout << "errsttotpr now is" << errsttotpr << endl;
  cout << endl;
  cout << "TOTAL PROMPT: $" << totpr << " \\pm " << sqrt(errsttotpr) << " \\pm " << sqrt(errsytotpr) << " (\\pm " << sqrt(errsttotpr+errsytotpr) << ")$" << endl;
  cout << "TOTAL NON-PROMPT: $" << totb << " \\pm " << sqrt(errsttotb) << " \\pm " << sqrt(errsytotb) << " (\\pm " << sqrt(errsttotb+errsytotb) << ")$" << endl;
  if (plotDellaBfrac) {

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TGraphAsymmErrors *o_al_y0 = new TGraphAsymmErrors(largeBins0, binAverages0, Bfracs0, binRMS0, binRMS0, Bfracerrs0, Bfracerrs0);
    TGraphAsymmErrors *o_al_y1 = new TGraphAsymmErrors(largeBins1, binAverages1, Bfracs1, binRMS1, binRMS1, Bfracerrs1, Bfracerrs1);
    TGraphAsymmErrors *o_al_y2 = new TGraphAsymmErrors(largeBins2, binAverages2, Bfracs2, binRMS2, binRMS2, Bfracerrs2, Bfracerrs2);

    o_al_y0->SetLineColor(kBlack);
    o_al_y0->SetMarkerStyle(21);
    o_al_y0->SetMarkerColor(kBlack);
    o_al_y1->SetLineColor(kBlue);
    o_al_y1->SetMarkerStyle(23);
    o_al_y1->SetMarkerColor(kBlue);
    o_al_y2->SetLineColor(kRed);
    o_al_y2->SetMarkerStyle(25);
    o_al_y2->SetMarkerColor(kRed);
    
    TMultiGraph *mg = new TMultiGraph("mg","mg");
    mg->Add(o_al_y0);
    mg->Add(o_al_y1);
    mg->Add(o_al_y2);
    
    TCanvas *d1 = new TCanvas("d1", "d1",600,600);
    // d1->cd(1)->SetLogy();
    d1->SetLeftMargin(0.15);
    mg->Draw("AP");
    o_al_y2->Draw("P");
    o_al_y1->Draw("P");
    o_al_y0->Draw("P");
    d1->Update();

    TH1F *hist=mg->GetHistogram();
    hist->SetXTitle("p_{T}^{J/#psi} (GeV/c)");
    hist->SetYTitle("b-fraction");
    
    TLatex * t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(13);
    t->SetTextFont(63);
    t->SetTextSizePixels(20);
    char lumi[100];
    sprintf(lumi, "CMS - #sqrt{s} = 7 TeV - L = 314 nb^{-1}");
    t->DrawLatex(0.48,0.4,lumi);
    
    TLegend *leg_a = new TLegend(0.66,0.16,0.92,0.35);
    o_al_y0->SetName("o_al_y0");
    o_al_y1->SetName("o_al_y1");
    o_al_y2->SetName("o_al_y2");
    leg_a->AddEntry("o_al_y0","|y|<1.2","p");
    leg_a->AddEntry("o_al_y1","1.2<|y|<1.6","p");
    leg_a->AddEntry("o_al_y2","1.6<|y|<2.4","p");
    leg_a->SetFillColor(0);
    leg_a->SetShadowColor(0);
    leg_a->SetLineColor(0);
    leg_a->Draw();
    
    d1->Print("allBfrac.pdf");
  }
  return;
}
