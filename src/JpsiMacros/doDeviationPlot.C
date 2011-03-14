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

void doDeviationPlot(string whichJpsi = "BJ", string whichPsip = "BP", string fileNameBase = "results", float minMax = 500., float minMaxPull = 10.) {

  static const unsigned int nbinspt1 = 6;  //6
  static const unsigned int nbinspt2 = 3;  //4
  static const unsigned int nbinspt3 = 6;  //8

  const float ptmin1 = 6.4;
  const float ptmin2 = 5.4;
  const float ptmin3 = 4.4;

  double ptbincenters1[nbinspt1] = {nbinspt1*20.};
  double ptbinerrors1[nbinspt1] = {nbinspt1*0.};
  double ptbincenters2[nbinspt2] = {nbinspt2*20.};
  double ptbinerrors2[nbinspt2] = {nbinspt2*0.};
  double ptbincenters3[nbinspt3] = {nbinspt3*20.};
  double ptbinerrors3[nbinspt3] = {nbinspt3*0.};

  double ycenters1[nbinspt1] = {nbinspt1*-999.};
  double yerrors1[nbinspt1] = {nbinspt1*1000000.};
  double pcenters1[nbinspt1] = {nbinspt1*-999.};
  double perrors1[nbinspt1] = {nbinspt1*1.};
  double ypcenters1[nbinspt1] = {nbinspt1*-999.};
  double yperrors1[nbinspt1] = {nbinspt1*1000000.};
  double ppcenters1[nbinspt1] = {nbinspt1*-999.};

  double ycenters2[nbinspt2] = {nbinspt2*-999.};
  double yerrors2[nbinspt2] = {nbinspt2*1000000.};
  double pcenters2[nbinspt2] = {nbinspt2*-999.};
  double perrors2[nbinspt2] = {nbinspt2*1.};
  double ypcenters2[nbinspt2] = {nbinspt2*-999.};
  double yperrors2[nbinspt2] = {nbinspt2*1000000.};
  double ppcenters2[nbinspt2] = {nbinspt2*-999.};

  double ycenters3[nbinspt3] = {nbinspt3*-999.};
  double yerrors3[nbinspt3] = {nbinspt3*1000000.};
  double pcenters3[nbinspt3] = {nbinspt3*-999.};
  double perrors3[nbinspt3] = {nbinspt3*1.};
  double ypcenters3[nbinspt3] = {nbinspt3*-999.};
  double yperrors3[nbinspt3] = {nbinspt3*1000000.};
  double ppcenters3[nbinspt3] = {nbinspt3*-999.};
  
  gROOT->ProcessLine(".! rm -f lista");
  char theCommand[50];
  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  
  ifstream lista("lista");
  char fileName[200];  
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  float trueMC, fitted, error; 
  int i = 0;         int j = 0;          int k = 0;
  string cutstring;
  float chi2i = 0.;  float chi2j = 0.;     float chi2k = 0.;
  float theMaximumPt1 = -1.;
  float theMinimumPt1 = 9000000.;
  float theMaximumPt2 = -1.;
  float theMinimumPt2 = 9000000.;
  float theMaximumPt3 = -1.;
  float theMinimumPt3 = 9000000.;

  // Inizia
  TH1F* pullJpsi = new TH1F("pullJpsi","Pulls for J/#psi",10,-6.,6.);
  TH1F* pullPsip = new TH1F("pullPsip","Pulls for #psi(2S)",10,-6.,6.);

  while (lista >> fileName) {

    cutstring = "results/" + fileNameBase + "_pT%f-%f_y%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      cout << "BIN " << ptmin << "-" << ptmax << " " << etamin << "-" << etamax << endl;

      ifstream theFile(fileName);
      if (fabs(etamin) < 0.001) {
	
	while (theFile >> theType >> trueMC >> fitted >> error) {

	  if (error < 0.000001) error = 0.000001;
	  if (ptmin > ptmin1 && !strcmp(theType,whichJpsi.c_str())) {

            if (ptmax > theMaximumPt1) theMaximumPt1 = ptmax;
            if (ptmin < theMinimumPt1) theMinimumPt1 = ptmin;
	    ptbincenters1[i] = (ptmax + ptmin)/2. ;
	    ptbinerrors1[i] = (ptmax - ptmin)/2. ; 
	    ycenters1[i] = fitted - trueMC;
	    yerrors1[i] = error;
            pcenters1[i] = ycenters1[i]/yerrors1[i];
            perrors1[i] = 1.0;
            chi2i += pow(ycenters1[i]/yerrors1[i],2);
	    pullJpsi->Fill(pcenters1[i]);
	    
	  } else if (ptmin > ptmin1 && !strcmp(theType,whichPsip.c_str())) {

	    ypcenters1[i] = fitted - trueMC;
	    yperrors1[i] = error;
            ppcenters1[i] = ypcenters1[i]/yperrors1[i];
            chi2i += pow(ycenters1[i]/yerrors1[i],2);
	    pullPsip->Fill(ppcenters1[i]);
	    i++;
	    
	  }
	}
      } else if (fabs(etamin-1.2) < 0.001) {
	
	while (theFile >> theType >> trueMC >> fitted >> error) {
	
	  if (error < 0.000001) error = 0.000001;
	  if (ptmin > ptmin2 && !strcmp(theType,whichJpsi.c_str())) {

            if (ptmax > theMaximumPt2) theMaximumPt2 = ptmax;
            if (ptmin < theMinimumPt2) theMinimumPt2 = ptmin;
	    ptbincenters2[j] = (ptmax + ptmin)/2. ;
	    ptbinerrors2[j] = (ptmax - ptmin)/2. ; 
	    ycenters2[j] = fitted - trueMC;
	    yerrors2[j] = error;
            pcenters2[j] = ycenters2[j]/yerrors2[j];
            perrors2[j] = 1.0;
            chi2j += pow(ycenters2[j]/yerrors2[j],2);
	    pullJpsi->Fill(pcenters2[j]);
	    
	  } else if (ptmin > ptmin2 && !strcmp(theType,whichPsip.c_str())) {

	    ypcenters2[j] = fitted - trueMC;
	    yperrors2[j] = error;
            ppcenters2[j] = ypcenters2[j]/yperrors2[j];
            chi2j += pow(ypcenters2[j]/yperrors2[j],2);
	    pullPsip->Fill(ppcenters2[j]);
	    j++;
	    
	  }
	}
      } else if (fabs(etamin-1.6) < 0.001) {
	
	while (theFile >> theType >> trueMC >> fitted >> error) {
	
	  if (error < 0.000001) error = 0.000001;
	  if (ptmin > ptmin3 && !strcmp(theType,whichJpsi.c_str())) {

            if (ptmax > theMaximumPt3) theMaximumPt3 = ptmax;
            if (ptmin < theMinimumPt3) theMinimumPt3 = ptmin;
	    ptbincenters3[k] = (ptmax + ptmin)/2. ;
	    ptbinerrors3[k] = (ptmax - ptmin)/2. ; 
	    ycenters3[k] = fitted - trueMC;
	    yerrors3[k] = error;
            pcenters3[k] = ycenters3[k]/yerrors3[k];
            perrors3[k] = 1.0;
            chi2k += pow(ycenters3[k]/yerrors3[k],2);
	    pullJpsi->Fill(pcenters3[k]);
	    
	  } else if (ptmin > ptmin3 && !strcmp(theType,whichPsip.c_str())) {

	    ypcenters3[k] = fitted - trueMC;
	    yperrors3[k] = error;
            ppcenters3[k] = ypcenters3[k]/yperrors3[k];
            chi2k += pow(ypcenters3[k]/yperrors3[k],2);
	    pullPsip->Fill(ppcenters3[k]);
	    k++;
	    
	  }
	}
      }

      theFile.close();
    }
  }

  lista.close();

  /* for (int aa = 0; aa < nbinspt; aa++) {
    cout << ptbincenters[aa] << " " << ptbinerrors[aa] << " " << ycenters1[aa] << " " << yerrors1[aa] << " " << ycenters2[aa] << " " << yerrors2[aa] << " " << endl;
    }*/
  // cout << theMinimumPt << " " << theMaximumPt << endl; 

  // Disegna residui
  TGraphErrors *gresi1 = new TGraphErrors(nbinspt1,ptbincenters1,ycenters1,ptbinerrors1,yerrors1);

  TCanvas c3("c3","c3",10,10,600,1000);
  c3.Divide(1,3);
  c3.cd(1);
  sprintf(theCommand,"Residual of fit results - barrel - %s",whichJpsi.c_str());
  gresi1->SetTitle(theCommand);
  gresi1->SetMarkerStyle(20);
  gresi1->SetMarkerColor(kRed);
  gresi1->SetMaximum(minMax);
  gresi1->SetMinimum(-minMax);
  gresi1->Draw("AP");

  TGraphErrors *presi1 = new TGraphErrors(nbinspt1,ptbincenters1,ypcenters1,ptbinerrors1,yperrors1);

  presi1->SetTitle(theCommand);
  presi1->SetMarkerStyle(21);
  presi1->SetMarkerColor(kBlue);
  presi1->Draw("PSAME");

  TLine aLine(theMinimumPt1,0.,theMaximumPt1,0.);
  aLine.SetLineStyle(kDashed);

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2i,nbinspt1));  
  TPaveLabel *pl1 = new TPaveLabel(0.2*(theMaximumPt1-theMinimumPt1),-0.9*minMax,
				   0.8*(theMaximumPt1-theMinimumPt1),-0.8*minMax,fileName);
  pl1->Draw("SAME");
  aLine.Draw("SAME");

  TGraphErrors *gresi2 = new TGraphErrors(nbinspt2,ptbincenters2,ycenters2,ptbinerrors2,yerrors2);

  c3.cd(2);
  sprintf(theCommand,"Residual of fit results - middle - %s",whichJpsi.c_str());
  gresi2->SetTitle(theCommand);
  gresi2->SetMarkerStyle(20);
  gresi2->SetMarkerColor(kRed);
  gresi2->SetMaximum(minMax);
  gresi2->SetMinimum(-minMax);
  gresi2->Draw("AP");

  TGraphErrors *presi2 = new TGraphErrors(nbinspt2,ptbincenters2,ypcenters2,ptbinerrors2,yperrors2);

  presi2->SetTitle(theCommand);
  presi2->SetMarkerStyle(21);
  presi2->SetMarkerColor(kBlue);
  presi2->Draw("PSAME");

  TLine aLine2(theMinimumPt2,0.,theMaximumPt2,0.);
  aLine2.SetLineStyle(kDashed);

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2j,nbinspt2));  
  TPaveLabel *pl2 = new TPaveLabel(0.2*(theMaximumPt2-theMinimumPt2),-0.9*minMax,
				   0.8*(theMaximumPt2-theMinimumPt2),-0.8*minMax,fileName);
  pl2->Draw("SAME");
  aLine2.Draw("SAME");

  TGraphErrors *gresi3 = new TGraphErrors(nbinspt3,ptbincenters3,ycenters3,ptbinerrors3,yerrors3);

  c3.cd(3);
  sprintf(theCommand,"Residual of fit results - forward - %s",whichJpsi.c_str());
  gresi3->SetTitle(theCommand);
  gresi3->SetMarkerStyle(20);
  gresi3->SetMarkerColor(kRed);
  gresi3->SetMaximum(minMax);
  gresi3->SetMinimum(-minMax);
  gresi3->Draw("AP");

  TGraphErrors *presi3 = new TGraphErrors(nbinspt3,ptbincenters3,ypcenters3,ptbinerrors3,yperrors3);

  presi3->SetTitle(theCommand);
  presi3->SetMarkerStyle(21);
  presi3->SetMarkerColor(kBlue);
  presi3->Draw("PSAME");

  TLine aLine3(theMinimumPt3,0.,theMaximumPt3,0.);
  aLine3.SetLineStyle(kDashed);

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2k,nbinspt3));  
  TPaveLabel *pl3 = new TPaveLabel(0.2*(theMaximumPt3-theMinimumPt3),-0.9*minMax,
				   0.8*(theMaximumPt3-theMinimumPt3),-0.8*minMax,fileName);
  pl3->Draw("SAME");
  aLine3.Draw("SAME");

  cutstring = "pictures/" + fileNameBase + "_" + whichJpsi + "resid.gif";
  c3.SaveAs(cutstring.c_str());

  // Disegna pull
  TGraphErrors *gpull1 = new TGraphErrors(nbinspt1,ptbincenters1,pcenters1,ptbinerrors1,perrors1);

  TCanvas c4("c4","c4",10,10,600,1000);
  c4.Divide(1,3);
  c4.cd(1);
  sprintf(theCommand,"Pull of fit results - barrel - %s",whichJpsi.c_str());
  gpull1->SetTitle(theCommand);
  gpull1->SetMarkerStyle(20);
  gpull1->SetMarkerColor(kRed);
  gpull1->SetMaximum(minMaxPull);
  gpull1->SetMinimum(-minMaxPull);
  gpull1->Draw("AP");

  TGraphErrors *ppull1 = new TGraphErrors(nbinspt1,ptbincenters1,ppcenters1,ptbinerrors1,perrors1);

  ppull1->SetTitle(theCommand);
  ppull1->SetMarkerStyle(21);
  ppull1->SetMarkerColor(kBlue);
  ppull1->Draw("PSAME");

  // pl1->Draw("SAME");
  aLine.Draw("SAME");

  TGraphErrors *gpull2 = new TGraphErrors(nbinspt2,ptbincenters2,pcenters2,ptbinerrors2,perrors2);

  c4.cd(2);
  sprintf(theCommand,"Pull of fit results - middle - %s",whichJpsi.c_str());
  gpull2->SetTitle(theCommand);
  gpull2->SetMarkerStyle(20);
  gpull2->SetMarkerColor(kRed);
  gpull2->SetMaximum(minMaxPull);
  gpull2->SetMinimum(-minMaxPull);
  gpull2->Draw("AP");

  TGraphErrors *ppull2 = new TGraphErrors(nbinspt2,ptbincenters2,ppcenters2,ptbinerrors2,perrors2);

  ppull2->SetTitle(theCommand);
  ppull2->SetMarkerStyle(21);
  ppull2->SetMarkerColor(kBlue);
  ppull2->Draw("PSAME");

  // pl2->Draw("SAME");
  aLine2.Draw("SAME");

  TGraphErrors *gpull3 = new TGraphErrors(nbinspt3,ptbincenters3,pcenters3,ptbinerrors3,perrors3);

  c4.cd(3);
  sprintf(theCommand,"Pull of fit results - forward - %s",whichJpsi.c_str());
  gpull3->SetTitle(theCommand);
  gpull3->SetMarkerStyle(20);
  gpull3->SetMarkerColor(kRed);
  gpull3->SetMaximum(minMaxPull);
  gpull3->SetMinimum(-minMaxPull);
  gpull3->Draw("AP");

  TGraphErrors *ppull3 = new TGraphErrors(nbinspt3,ptbincenters3,ppcenters3,ptbinerrors3,perrors3);

  ppull3->SetTitle(theCommand);
  ppull3->SetMarkerStyle(21);
  ppull3->SetMarkerColor(kBlue);
  ppull3->Draw("PSAME"); 

  // pl3->Draw("SAME");
  aLine3.Draw("SAME");

  cutstring = "pictures/" + fileNameBase + "_" + whichJpsi + "pull.gif";
  c4.SaveAs(cutstring.c_str());

  // Istogramma pull
  TCanvas c5("c5","c5",10,10,600,400);
  c5.cd();
  pullPsip->SetLineColor(kRed);
  pullPsip->SetMarkerColor(kRed);
  pullPsip->Fit("gaus","","e1");
  pullJpsi->Fit("gaus","","e1same");

  cutstring = "pictures/" + fileNameBase + "_" + whichJpsi + "pullHist.gif";
  c5.SaveAs(cutstring.c_str());

  return;

  
}
