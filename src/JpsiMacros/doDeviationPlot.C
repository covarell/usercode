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

void doDeviationPlot(string whichType = "GG", string fileNameBase = "results", float minMax = 500.) {

  static const unsigned int nbinspt = 7;

  double ptbincenters[nbinspt] = {nbinspt*20.};
  double ptbinerrors[nbinspt] = {nbinspt*0.};
  double ycenters1[nbinspt] = {nbinspt*-999.};
  double yerrors1[nbinspt] = {nbinspt*1000000.};
  double ycenters2[nbinspt] = {nbinspt*-999.};
  double yerrors2[nbinspt] = {nbinspt*1000000.};
  double ycenters3[nbinspt] = {nbinspt*-999.};
  double yerrors3[nbinspt] = {nbinspt*1.};
  double ycenters4[nbinspt] = {nbinspt*-999.};
  
  gROOT->ProcessLine(".! rm -f lista");
  char theCommand[50];
  sprintf(theCommand,".! ls results/%s*.txt > lista",fileNameBase.c_str());
  gROOT->ProcessLine(theCommand);
  
  ifstream lista("lista");
  char fileName[200];  
  float ptmin, ptmax, etamin, etamax; 
  char theType[2];
  float trueMC, fitted, error; 
  int i = 0;         int j = 0;
  string cutstring;
  float chi2i = 0.;  float chi2j = 0.;
  float theMaximumPt = -1.;
  float theMinimumPt = 9000000.;

  // Inizia

  while (lista >> fileName) {

    cutstring = "results/" + fileNameBase + "_pT%f-%f_eta%f-%f.txt"; 
    if (strstr(fileName,fileNameBase.c_str()) && sscanf(fileName, cutstring.c_str(), &ptmin, &ptmax, &etamin, &etamax) != 0) {

      ifstream theFile(fileName);
      if (etamin == 0.0) {
	
	while (theFile >> theType >> trueMC >> fitted >> error) {
	
	  if (!strcmp(theType,whichType.c_str())) {

            if (ptmax > theMaximumPt) theMaximumPt = ptmax;
            if (ptmin < theMinimumPt) theMinimumPt = ptmin;
	    ptbincenters[i] = (ptmax + ptmin)/2. ;
	    ptbinerrors[i] = (ptmax - ptmin)/2. ; 
	    ycenters1[i] = fitted - trueMC;
	    yerrors1[i] = error;
            ycenters3[i] = ycenters1[i]/yerrors1[i];
            yerrors3[i] = 1.0;
            chi2i += pow(ycenters1[i]/yerrors1[i],2);
	    i++;
	    
	  }
	}  
      } else if (etamax == 2.5) {
	while (theFile >> theType >> trueMC >> fitted >> error) {

	  if (!strcmp(theType,whichType.c_str())) {
            if (ptmax > theMaximumPt) theMaximumPt = ptmax;
            if (ptmin < theMinimumPt) theMinimumPt = ptmin;
	    ycenters2[j] = fitted - trueMC;
	    yerrors2[j] = error;
            ycenters4[j] = ycenters2[j]/yerrors2[j];
            yerrors3[j] = 1.0;
            chi2j += pow(ycenters2[j]/yerrors2[j],2);
	    j++;
	  }
	} 
      }

      theFile.close();
    }
  }

  lista.close();

  for (int aa = 0; aa < nbinspt; aa++) {
    cout << ptbincenters[aa] << " " << ptbinerrors[aa] << " " << ycenters1[aa] << " " << yerrors1[aa] << " " << ycenters2[aa] << " " << yerrors2[aa] << " " << endl;
  }
  // cout << theMinimumPt << " " << theMaximumPt << endl; 

  TGraphErrors *gpull1 = new TGraphErrors(nbinspt,ptbincenters,ycenters1,ptbinerrors,yerrors1);

  TCanvas c3;
  c3.cd();
  sprintf(theCommand,"Residual of fit results - barrel - %s",whichType.c_str());
  gpull1->SetTitle(theCommand);
  gpull1->SetMarkerStyle(20);
  gpull1->SetMarkerColor(kRed);
  gpull1->SetMaximum(minMax);
  gpull1->SetMinimum(-minMax);  if (!strcmp("RE",whichType.c_str())) gpull1->SetMinimum(2.);
  gpull1->Draw("AP");

  TLine aLine(theMinimumPt,0.,theMaximumPt,0.);
  aLine.SetLineStyle(kDashed);

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2i,nbinspt));  
  TPaveLabel *pl1 = new TPaveLabel(0.2*(theMaximumPt-theMinimumPt),-0.9*minMax,
				   0.8*(theMaximumPt-theMinimumPt),-0.8*minMax,fileName);
  if (strcmp("RE",whichType.c_str())) {
    pl1->Draw("SAME");
    aLine.Draw("SAME");
  }

  cutstring = "pictures/" + fileNameBase + "_" + whichType + "resid_barrel.gif";
  c3.SaveAs(cutstring.c_str());

  TGraphErrors *gpull2 = new TGraphErrors(nbinspt,ptbincenters,ycenters2,ptbinerrors,yerrors2);

  TCanvas c4;
  c4.cd();
  sprintf(theCommand,"Residual of fit results - endcap - %s",whichType.c_str());
  gpull2->SetTitle(theCommand);
  gpull2->SetMarkerStyle(20);
  gpull2->SetMarkerColor(kBlue);
  gpull2->SetMaximum(minMax);
  gpull2->SetMinimum(-minMax);   if (!strcmp("RE",whichType.c_str())) gpull2->SetMinimum(2.);
  gpull2->Draw("AP");

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2j,nbinspt));  
  TPaveLabel *pl2 = new TPaveLabel(0.2*(theMaximumPt-theMinimumPt),-0.9*minMax,
				   0.8*(theMaximumPt-theMinimumPt),-0.8*minMax,fileName);
  if (strcmp("RE",whichType.c_str())) { 
    pl2->Draw("SAME");
    aLine.Draw("SAME");
  }

  cutstring = "pictures/" + fileNameBase + "_" + whichType + "resid_endcap.gif";
  c4.SaveAs(cutstring.c_str());

  TGraphErrors *gpull3 = new TGraphErrors(nbinspt,ptbincenters,ycenters3,ptbinerrors,yerrors3);

  c3.cd();
  sprintf(theCommand,"Pull of fit results - barrel - %s",whichType.c_str());
  gpull3->SetTitle(theCommand);
  gpull3->SetMarkerStyle(20);
  gpull3->SetMarkerColor(kRed);
  gpull3->SetMaximum(5.);
  gpull3->SetMinimum(-5.);
  gpull3->Draw("AP");

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2i,nbinspt));  
  TPaveLabel *pl1 = new TPaveLabel(0.2*(theMaximumPt-theMinimumPt),-4.5,
				   0.8*(theMaximumPt-theMinimumPt),-3.9,fileName);

  pl1->Draw("SAME");
  aLine.Draw("SAME");

  cutstring = "pictures/" + fileNameBase + "_" + whichType + "pull_barrel.gif";
  c3.SaveAs(cutstring.c_str());

  TGraphErrors *gpull4 = new TGraphErrors(nbinspt,ptbincenters,ycenters4,ptbinerrors,yerrors3);

  c4.cd();
  sprintf(theCommand,"Residual of fit results - endcap - %s",whichType.c_str());
  gpull4->SetTitle(theCommand);
  gpull4->SetMarkerStyle(20);
  gpull4->SetMarkerColor(kBlue);
  gpull4->SetMaximum(5.);
  gpull4->SetMinimum(-5.);
  gpull4->Draw("AP");

  sprintf(fileName,"Pearson's test prob. = %f",TMath::Prob(chi2j,nbinspt));  
  TPaveLabel *pl2 = new TPaveLabel(0.2*(theMaximumPt-theMinimumPt),-4.5,
				   0.8*(theMaximumPt-theMinimumPt),-3.9,fileName);

  pl2->Draw("SAME");
  aLine.Draw("SAME");

  cutstring = "pictures/" + fileNameBase + "_" + whichType + "pull_endcap.gif";
  c4.SaveAs(cutstring.c_str());
        
  return;
}
