//===================================================================
//H.K. Woehri and C.Lourenco [25. July 2010]
//LHC data presented at the ICHEP conference, 22-28 July 2010, Paris
//CMS data from the paper, based on 314 nb-1
//===================================================================
// #include "lhcSettings.inc"

Int_t const colourLHCb = 2;
Int_t const colourALICE = kMagenta-3;
Int_t const colourCMS[2] = {kCyan-3, kBlue};
Int_t const colourATLAS[3] = {3,kGreen-1,3};
//
Int_t const markerCMS[2] = {21,20}; 
Int_t const markerLHCb = 29;
Int_t const markerALICE = 25;
Int_t const markerATLAS[3] = {22,23,22};

Int_t const kNbRap = 3;
Double_t rapRange[kNbRap+1] = {0., 1.2, 1.6, 2.4};
Int_t const kNbpTMid = 3;
Int_t const kNbpTInter = 4;
Int_t const kNbpTFWD = 8;
TGraphAsymmErrors *gB_mid, *gB_inter, *gB_FWD;
TGraphErrors *gB_mid_syst, *gB_inter_syst, *gB_FWD_syst;
TGraphErrors *gB_CDF_syst, *gB_CDF;
TGraphErrors *gB_LHCb, *gB_LHCb_differential, *gB_LHCb_differential_syst;
TGraphErrors *gB_ATLAS, *gB_ATLAS_syst;

void LoadCMS();
void LoadLHCb();
void LoadATLAS();
void LoadCDF();
void PlotData();

//======================
void bFraction(){

  LoadCMS();
  LoadCDF();
  LoadLHCb();
  LoadATLAS();

  PlotCMS();
}
//======================
void PlotCMS(){

  Double_t dummyX[1] = {0.}, dummyY[1] = {0.};
  TGraph *gDummy = new TGraph(1, dummyX, dummyY);

  Double_t fontSize = 0.04;
  Double_t maxX = 19.9;
  //===================================
  //CMS alone
  //===================================
  TCanvas *c1 = new TCanvas("c1cms", "c1");
  TH1F *hFrame1 = gPad->DrawFrame(0., 0., maxX, 0.55);
  hFrame1->SetXTitle("p_{T}^{J/#psi} [GeV/c]");
  hFrame1->SetYTitle("Fraction of J/#psi from b hadrons");

  //redefine the colours:
  gB_mid_syst->SetLineColor(3);
  gB_inter_syst->SetLineColor(4);
  gB_FWD_syst->SetLineColor(2);
  gB_FWD_syst->SetMarkerColor(2);

  gB_mid->SetLineColor(3);   gB_mid->SetMarkerColor(3);
  // gB_mid->SetMarkerSize(2*gB_mid->GetMarkerSize());
  gB_inter->SetLineColor(4); gB_inter->SetMarkerColor(4); gB_FWD->SetMarkerStyle(22);  
  gB_FWD->SetLineColor(2);   gB_FWD->SetMarkerColor(2);  
  //
  gB_mid_syst->Draw("e2");
  gB_inter_syst->Draw("e2");
  gB_FWD_syst->Draw("e2");

  gB_mid->Draw("p");
  gB_inter->Draw("p");
  gB_FWD->Draw("p");

  // TLatex *tex1 = new TLatex(0.6887088,0.34, "#sqrt{s} = 7 TeV");
  // tex1->SetTextSize(fontSize+0.01);
  // tex1->Draw();

  TLegend *leg1 = new TLegend(0.1521839,0.720339,0.3505747,0.9427966, "CMS   #sqrt{s} = 7 TeV   314 nb^{-1}");
  leg1->AddEntry(gB_FWD, "1.6 < |y| < 2.4", "p");
  leg1->AddEntry(gB_inter, "1.2 < |y| < 1.6", "p");
  leg1->AddEntry(gB_mid, "|y| < 1.2", "p");
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(fontSize+0.005); leg1->Draw();

  c1->Print("Figures/bFraction_CMS_314nb.eps");
  c1->Print("Figures/bFraction_CMS_314nb.pdf");
  c1->Print("Figures/bFraction_CMS_314nb.gif");

  //===================================
  //CMS with CDF
  //===================================
  TCanvas *c2 = new TCanvas("c2cms", "c2");
  TH1F *hFrame2 = gPad->DrawFrame(0., 0., maxX, 0.55);
  hFrame2->SetXTitle("p_{T}^{J/#psi} [GeV/c]");
  hFrame2->SetYTitle("Fraction of J/#psi from b hadrons");

  //redefine the colours:
  gB_mid_syst->SetLineColor(3);
  gB_inter_syst->SetLineColor(4);
  gB_FWD_syst->SetLineColor(2);

  gB_mid->SetLineColor(3);   gB_mid->SetMarkerColor(3);
  gB_inter->SetLineColor(4); gB_inter->SetMarkerColor(4); gB_FWD->SetMarkerStyle(22);  
  gB_FWD->SetLineColor(2);   gB_FWD->SetMarkerColor(2);  
  //
  gB_mid_syst->Draw("p");
  gB_inter_syst->Draw("p");
  gB_FWD_syst->Draw("p");

  gB_CDF_syst->Draw("p");

  gB_mid->Draw("p");
  gB_inter->Draw("p");
  gB_FWD->Draw("p");

  gB_CDF->Draw("p");

  // TLatex *tex2 = new TLatex(0.6887088,0.34, "#sqrt{s} = 7 TeV");
  // tex2->SetTextSize(fontSize+0.01);
  // tex2->Draw();

  TLegend *leg2 = new TLegend(0.1821839,0.710339,0.3505747,0.9327966, "CMS   #sqrt{s} = 7 TeV   314 nb^{-1}");
  leg2->AddEntry(gB_FWD, "1.6 < |y| < 2.4", "p");
  leg2->AddEntry(gB_inter, "1.2 < |y| < 1.6", "p");
  leg2->AddEntry(gB_mid, "|y| < 1.2", "p");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(fontSize+0.005); leg2->Draw();

  TLatex *tex2 = new TLatex(8.00517,0.02560343, "PRD 71 (2005) 032001");
  tex2->SetTextSize(fontSize);
  tex2->Draw();

  TLegend *leg2 = new TLegend(0.3531609,0.2055085,0.762931,0.2966102);
  leg2->AddEntry(gB_CDF, "CDF   #sqrt{s} = 1.96 TeV   |y| < 0.6", "p");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(fontSize); leg2->Draw();

  c2->Print("Figures/bFraction_CMS_314nb_CDF.eps");
  c2->Print("Figures/bFraction_CMS_314nb_CDF.pdf");
  c2->Print("Figures/bFraction_CMS_314nb_CDF.gif");
}

//======================
void LoadCMS(){

  //================================
  //BPH-10-002, 314 nb-1
  //================================
  //mid-rapidity, |y| < 1.2
  Double_t bFraction_mid[kNbpTMid] = {0.179, 0.257, 0.395};
  Double_t bFrac_Stat_mid[kNbpTMid] = {0.024, 0.015, 0.018};
  Double_t bFrac_Syst_mid[kNbpTMid] = {0.015, 0.014, 0.005};
  Double_t bFrac_tot_mid[kNbpTMid];
  Double_t pTMid[kNbpTMid] = {5.46, 8.14, 13.5};
  Double_t pTMid_Err[kNbpTMid] = {0., 0., 0.};
  Double_t pTRange_mid[kNbpTMid+1] = {4.5,6.5,10.,30.};
  Double_t pTL_mid[kNbpTMid], pTR_mid[kNbpTMid];
  Double_t errPTSyst_mid[kNbpTMid];
  for(int ipT = 0; ipT < kNbpTMid; ipT++){
    // pTL_mid[ipT] = 0.;
    // pTR_mid[ipT] = 0.;
    pTL_mid[ipT] = pTMid_Err[ipT];
    pTR_mid[ipT] = pTMid_Err[ipT];
    errPTSyst_mid[ipT] = 0.;
    bFrac_tot_mid[ipT] = sqrt(pow(bFrac_Stat_mid[ipT],2) + pow(bFrac_Syst_mid[ipT],2));
  }
  gB_mid = new TGraphAsymmErrors(kNbpTMid, pTMid, bFraction_mid, pTL_mid, pTR_mid, bFrac_tot_mid,bFrac_tot_mid);
  // gB_mid->SetMarkerStyle(markerCMS[0]);
  gB_mid->SetMarkerStyle(markerCMS[0]);
  gB_mid->SetMarkerColor(colourCMS[0]);
  gB_mid->SetLineColor(colourCMS[0]);
  gB_mid->SetMarkerSize(1.1);
  gB_mid_syst = new TGraphErrors(kNbpTMid, pTMid, bFraction_mid, errPTSyst_mid, bFrac_Stat_mid);
  gB_mid_syst->SetLineColor(colourCMS[0]);
  gB_mid_syst->SetFillStyle(0);

  //intermediate-rapidity, 1.2 < |y| < 1.6
  Double_t bFraction_inter[kNbpTInter] = {0.146, 0.180, 0.203, 0.360};
  Double_t bFrac_Stat_inter[kNbpTInter] = {0.021, 0.017, 0.017, 0.031};
  Double_t bFrac_Syst_inter[kNbpTInter] = {0.028, 0.019, 0.014, 0.016};
  Double_t bFrac_tot_inter[kNbpTInter];
  Double_t pTInter[kNbpTInter] = {3.27, 5.48, 7.89, 12.96};
  Double_t pTInter_Err[kNbpTInter] = {0., 0., 0., 0.};
  Double_t pTRange_inter[kNbpTInter+1] = {2.0,4.5,6.5,10.,30.};
  Double_t pTL_inter[kNbpTInter], pTR_inter[kNbpTInter];
  Double_t errPTSyst_inter[kNbpTInter];
  for(int ipT = 0; ipT < kNbpTInter; ipT++){
    // pTL_inter[ipT] = 0.;
    // pTR_inter[ipT] = 0.;
    pTL_inter[ipT] = pTInter_Err[ipT];
    pTR_inter[ipT] = pTInter_Err[ipT];
    errPTSyst_inter[ipT] = 0.;
    bFrac_tot_inter[ipT] = sqrt(pow(bFrac_Stat_inter[ipT],2) + pow(bFrac_Syst_inter[ipT],2));
  }
  gB_inter = new TGraphAsymmErrors(kNbpTInter, pTInter, bFraction_inter, pTL_inter, pTR_inter, bFrac_tot_inter, bFrac_tot_inter);
  gB_inter->SetMarkerStyle(markerCMS[1]);
  gB_inter->SetMarkerColor(colourCMS[1]);
  gB_inter->SetLineColor(colourCMS[1]);
  gB_inter->SetMarkerSize(1.1);
  gB_inter_syst = new TGraphErrors(kNbpTInter, pTInter, bFraction_inter, errPTSyst_inter, bFrac_Stat_inter);
  gB_inter_syst->SetLineColor(colourCMS[1]);
  gB_inter_syst->SetFillStyle(0);

  //forward-rapidity, 1.6 < |y| < 2.4
  Double_t bFraction_FWD[kNbpTFWD] = {0.057, 0.087, 0.113, 0.139, 0.160, 0.177, 0.235, 0.374};
  Double_t bFrac_Stat_FWD[kNbpTFWD] = {0.021, 0.014, 0.013, 0.014, 0.014, 0.012, 0.016, 0.031};
  Double_t bFrac_Syst_FWD[kNbpTFWD] = {0.042, 0.022, 0.020, 0.010, 0.013, 0.012, 0.012, 0.008};
  Double_t bFrac_tot_FWD[kNbpTFWD];
  Double_t pT_FWD[kNbpTFWD] = {0.79, 1.60, 2.35, 3.10, 3.96, 5.35, 7.86, 13.11};
  Double_t pT_FWD_Err[kNbpTFWD] = {0., 0., 0., 0., 0., 0., 0., 0.};
  Double_t pTRange_FWD[kNbpTFWD+1] = {0., 1.25, 2., 2.75, 3.5, 4.5, 6.5, 10., 30.};
  Double_t pTL_FWD[kNbpTFWD], pTR_FWD[kNbpTFWD];
  Double_t errPTSyst_FWD[kNbpTFWD];
  for(int ipT = 0; ipT < kNbpTFWD; ipT++){
    // pTL_FWD[ipT] = 0.;
    // pTR_FWD[ipT] = 0.;
    pTL_FWD[ipT] = pT_FWD_Err[ipT];
    pTR_FWD[ipT] = pT_FWD_Err[ipT];
    errPTSyst_FWD[ipT] = 0.;
    bFrac_tot_FWD[ipT] = sqrt(pow(bFrac_Stat_FWD[ipT],2) + pow(bFrac_Syst_FWD[ipT],2));
  }
  gB_FWD = new TGraphAsymmErrors(kNbpTFWD, pT_FWD, bFraction_FWD, pTL_FWD, pTR_FWD, bFrac_tot_FWD, bFrac_tot_FWD);
  gB_FWD->SetMarkerStyle(markerCMS[2]);
  gB_FWD->SetMarkerColor(colourCMS[2]);
  gB_FWD->SetLineColor(colourCMS[2]);
  gB_FWD->SetMarkerSize(1.1);
  gB_FWD_syst = new TGraphErrors(kNbpTFWD, pT_FWD, bFraction_FWD, errPTSyst_FWD, bFrac_Stat_FWD);
  gB_FWD_syst->SetLineColor(colourCMS[2]);
  gB_FWD_syst->SetMarkerStyle(markerCMS[2]);
  gB_FWD_syst->SetMarkerColor(colourCMS[2]);  
  gB_FWD_syst->SetFillStyle(0);
  
}

//======================
void LoadCDF(){

  //================================
  //Phys. Rev. D71 (2005) 032001
  //================================
  Int_t const kNbPoints = 26;
  Double_t bFraction[kNbPoints] = {0.094, 0.092, 0.085, 0.100,
				   0.091, 0.101, 0.099, 0.109,
				   0.112, 0.113, 0.133, 0.116,
				   0.126, 0.131, 0.147, 0.141,
				   0.156, 0.169, 0.182, 0.208,
				   0.227, 0.250, 0.279, 0.337,
				   0.397, 0.464};
  Double_t bFrac_Stat[kNbPoints] = {0.010, 0.006, 0.006, 0.005,
				   0.005, 0.005, 0.005, 0.005,
				   0.005, 0.005, 0.005, 0.005,
				   0.006, 0.006, 0.007, 0.005,
				   0.006, 0.007, 0.007, 0.006,
				   0.009, 0.011, 0.012, 0.019,
				   0.025, 0.045};
  Double_t bFrac_Syst[kNbPoints] = {0.012, 0.010, 0.009, 0.011,
				   0.010, 0.009, 0.008, 0.007,
				   0.008, 0.007, 0.007, 0.007,
				   0.007, 0.007, 0.008, 0.006,
				   0.007, 0.007, 0.008, 0.009,
				   0.007, 0.008, 0.008, 0.009,
				   0.009, 0.014};
  Double_t pTMid[kNbPoints] = {1.38, 1.63, 1.87, 2.13, 2.38, 2.62, 
			      2.87, 3.12, 3.38, 3.62, 3.87, 4.12, 
			      4.38, 4.62, 4.88, 5.24, 5.75, 6.24, 
			      6.74, 7.45, 8.46, 9.46, 10.8, 12.8, 
			      15.2, 18.3};
  Double_t bFrac_tot[kNbPoints];
  Double_t errPT[kNbPoints];
  Double_t errPTSyst[kNbPoints];
  for(int ipT = 0; ipT < kNbPoints; ipT++){
    errPT[ipT] = 0.;
    errPTSyst[ipT] = 0.;
    bFrac_tot[ipT] = sqrt(pow(bFrac_Stat[ipT],2) + pow(bFrac_Syst[ipT],2));
  }
  gB_CDF = new TGraphErrors(kNbPoints, pTMid, bFraction, errPT, bFrac_tot);
  gB_CDF->SetMarkerStyle(20);
  gB_CDF->SetMarkerSize(0.7);
  gB_CDF_syst = new TGraphErrors(kNbPoints, pTMid, bFraction, errPTSyst, bFrac_Stat);
  gB_CDF_syst->SetFillColor(11);
  gB_CDF_syst->SetMarkerSize(0.1); 
}

//==========================
void LoadLHCb(){

  //==============================
  //values from ICHEP presentation
  //(pT integrated)
  //==============================
  Int_t const kNbPoints = 10;
  Double_t dNdpT[kNbPoints] = {427., 823., 687., 398., 259., 
			    163., 74., 34., 23., 10.};
  Double_t dNdpTErr[kNbPoints] = {31., 40., 36., 24., 18.,
				  13., 9., 6., 5., 3.};
  Double_t pT[kNbPoints];
  for(int iP = 0; iP < kNbPoints; iP++)
    pT[iP] = 0.5 + iP;
  
  Double_t weightedAv = 0., weight = 0.;
  for(int iP = 0; iP < kNbPoints; iP++){

    weightedAv += pT[iP] * dNdpT[iP];
    weight += dNdpT[iP];
  }
  weightedAv /= weight;
  printf("average pT: %1.3f GeV/c\n", weightedAv);

  //value, averaged over all pT:
  Double_t fracB[1] = {0.111};
  Double_t fracB_err[1] = {0.008};
  Double_t pTMean[1] = {weightedAv};
  Double_t errPT[1] = {0.};

  gB_LHCb = new TGraphErrors(1, pTMean, fracB, errPT, fracB_err);
  gB_LHCb->SetMarkerStyle(markerLHCb);
  gB_LHCb->SetMarkerSize(1.4);
  gB_LHCb->SetLineWidth(2);
  gB_LHCb->SetMarkerColor(colourLHCb);
  gB_LHCb->SetLineColor(colourLHCb);

  //==============================
  //values sent by Patrick Robbe
  //mail from 24 July 2010
  //--> Conf. report ???
  //==============================
  Double_t pTAv[4] = {1.23, 2.87, 4.86, 7.24};
  Double_t fracB_LHCb[4] = {0.091, 0.098, 0.150, 0.240};
  Double_t fracB_LHCb_err[4] = {0.013, 0.013, 0.022, 0.046};
  Double_t fracB_LHCb_syst[4] = {0.005, 0.005, 0.008, 0.012};
  Double_t errPT_LHCb[4] = {0., 0., 0., 0.};
  Double_t errPT_LHCb_syst[4];
  for(int iP = 0; iP < 4; iP++)
    errPT_LHCb_syst[iP] = 0.1;

  gB_LHCb_differential = new TGraphErrors(4, pTAv, fracB_LHCb, errPT_LHCb, fracB_LHCb_err);
  gB_LHCb_differential->SetMarkerStyle(markerLHCb);
  gB_LHCb_differential->SetMarkerSize(1.4);
  gB_LHCb_differential->SetLineWidth(2);
  gB_LHCb_differential->SetMarkerColor(colourLHCb);
  gB_LHCb_differential->SetLineColor(colourLHCb);

  gB_LHCb_differential_syst = new TGraphErrors(4, pTAv, fracB_LHCb, errPT_LHCb_syst, fracB_LHCb_syst);
  gB_LHCb_differential_syst->SetLineColor(colourLHCb);
  gB_LHCb_differential_syst->SetFillStyle(0);


}
//==========================
void LoadATLAS(){

  //=================================
  //data presented at ICHEP
  //by Andrew Nelson
  //derived from ratios shown on p. 21
  //==================================

  Double_t ratioBovP[5] = {0.22, 0.12, 0.24, 0.25, 0.60};
  Double_t ratioBovP_stat[5] = {0.09, 0.05, 0.05, 0.08, 0.15};
  Double_t ratioBovP_syst[5] = {0.07, 0.06, 0.05, 0.07, 0.10};
  Double_t pT[5] = {2.5, 5.0, 7.0, 9.0, 12.5};
  Double_t errPT[5] = {0., 0., 0., 0., 0.};

  Double_t fracB[5], fracB_err[5], fracB_syst[5];
  Double_t errPTSyst[5];
  for(int iP = 0; iP < 5; iP++){
    fracB[iP] = ratioBovP[iP] / (ratioBovP[iP] + 1.);
    fracB_err[iP] = ratioBovP_stat[iP] / pow(ratioBovP[iP]+1.,2);
    fracB_syst[iP] = ratioBovP_syst[iP] / pow(ratioBovP[iP]+1.,2);
    errPTSyst[iP] = 0.1;
  }

  gB_ATLAS = new TGraphErrors(5, pT, fracB, errPT, fracB_err);
  gB_ATLAS->SetMarkerStyle(markerATLAS[0]);
  gB_ATLAS->SetMarkerSize(1.2);
//   gB_ATLAS->SetLineWidth(2);
  gB_ATLAS->SetMarkerColor(colourATLAS[0]);
  gB_ATLAS->SetLineColor(colourATLAS[0]);

  gB_ATLAS_syst = new TGraphErrors(5, pT, fracB, errPTSyst, fracB_syst);
  gB_ATLAS_syst->SetLineColor(colourATLAS[0]);
  gB_ATLAS_syst->SetFillStyle(0);

}
