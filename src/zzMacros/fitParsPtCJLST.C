// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2F.h>

static const unsigned int Nmasspoints = 9;
string masspointsS[Nmasspoints] = {"115","120","125","130","140","200","400","700","1000"};

// static const unsigned int Nmasspoints = 7;
// string masspointsS[Nmasspoints] = {"120","130","140","200","400","700","1000"};
string nameSample[8] = {"gg","vbf","zz","zx","ggzz","wh","zh","tth"};
TString systSources[8][5];
bool isProfiledSource[8][5];  // is a source is not "profiled" with mass
                              // take relative difference at 125 and apply
                              // everywhere

Double_t zzFunc(Double_t *x, Double_t *par) {
  Double_t fitval;
  if (x[0] <= 182.) fitval = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  else fitval = par[3]*x[0] + par[4]*x[0]*x[0] - par[3]*182. - par[4]*33124. + par[0] + par[1]*182. + par[2]*33124.;
  return fitval;
}

void fitAndWrite(Double_t *x, Double_t *y, Double_t *xerr, Double_t *yerr, string parName, TF1* funcName, ofstream *allParams) {

  int nPointsRemoved = 0;
  TGraphErrors *yfit = new TGraphErrors(Nmasspoints,x,y,xerr,yerr);

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    if (fabs(y[i]) > 999999.) {
      cout << parName << ": Point of mass " << masspointsS[i] << " is not there, removed" << endl; 
      yfit->RemovePoint(i-nPointsRemoved);
      nPointsRemoved++;
    }
  }
  
  yfit->SetTitle(parName.c_str());
  yfit->SetMarkerStyle(20);
  yfit->SetMarkerColor(kRed);
  yfit->Draw("AP"); 
  // gPad->SetLogx(); 
  if (yerr[1] > 0.0000001 || yerr[Nmasspoints-1] > 0.0000001) {
    yfit->Fit(funcName,"","PSAME"); 
    for (unsigned int j = 0; j < 5; j++) {
      *allParams << parName << " [" << j << "] = " << funcName->GetParameter(j) << "\n";
    }
  } else {
    TF1 *mypol0 = new TF1("mypol0","pol0");
    yfit->Fit("mypol0","","PSAME"); 
    *allParams << parName << " [0] = " << mypol0->GetParameter(0) << "\n";
  }
}

void fitParsPtCJLST(int LHCsqrts = 7, int whichtype = 1) {

// whichtype
// 0 - gg Signal
// 1 - VBF Signal
// 2 - ZZ
// 3 - ZX
// 4 - ggZZ
// 5 - WH
// 6 - ZH
// 7 - ttH  

  cout << "FIRST RUN source changeCinZero.csh!" << endl;

  systSources[0][0] = "Resummation";
  systSources[0][1] = "TopMass";
  systSources[0][2] = "Mela";
  
  systSources[1][0] = "PDF-VBF";
  systSources[1][1] = "scale-VBF";
  systSources[1][2] = "Mela";
  
  systSources[2][0] = "NLOLO_WH";
  
  systSources[3][0] = "NLOLO_ZH";

  systSources[5][0] = "SingleZ";
  systSources[5][1] = "PDF-ZZ";
  systSources[5][2] = "scale-ZZ";
  systSources[5][3] = "Mela";

  isProfiledSource[0][0] = false;
  isProfiledSource[0][1] = false;
  isProfiledSource[0][2] = false;
  
  isProfiledSource[1][0] = false;
  isProfiledSource[1][1] = false;
  isProfiledSource[1][2] = false;
  
  isProfiledSource[2][0] = false;
  
  isProfiledSource[3][0] = false;

  isProfiledSource[5][0] = false;
  isProfiledSource[5][1] = false;
  isProfiledSource[5][2] = false;
  isProfiledSource[5][3] = false;

  double masspoints[Nmasspoints] = {Nmasspoints*0.};
  double masspointserr[Nmasspoints] = {Nmasspoints*0.};

  double bb[Nmasspoints] = {Nmasspoints*1000000.}; 
  double T[Nmasspoints] = {Nmasspoints*1000000.}; 
  double n[Nmasspoints] = {Nmasspoints*1000000.}; 
  double bbdue[Nmasspoints] = {Nmasspoints*1000000.}; 
  double ndue[Nmasspoints] = {Nmasspoints*1000000.}; 
  double m[Nmasspoints] = {Nmasspoints*1000000.}; 
  double fexp[Nmasspoints] = {Nmasspoints*1000000.}; 

  double bberr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double Terr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double nerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double bbdueerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double ndueerr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double merr[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double fexperr[Nmasspoints] = {Nmasspoints*0.000000001};

  double bberrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double Terrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double nerrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double bbdueerrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double ndueerrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double merrup[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double fexperrup[Nmasspoints] = {Nmasspoints*0.000000001};

  double bberrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double Terrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double nerrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double bbdueerrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double ndueerrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double merrdown[Nmasspoints] = {Nmasspoints*0.000000001}; 
  double fexperrdown[Nmasspoints] = {Nmasspoints*0.000000001};

  double bbup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double Tup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double nup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double bbdueup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double ndueup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double mup[Nmasspoints] = {Nmasspoints*1000000.}; 
  double fexpup[Nmasspoints] = {Nmasspoints*1000000.}; 

  double bbdown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double Tdown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double ndown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double bbduedown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double nduedown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double mdown[Nmasspoints] = {Nmasspoints*1000000.}; 
  double fexpdown[Nmasspoints] = {Nmasspoints*1000000.}; 

  ofstream *allParams;
  ofstream *allParamsUp;
  ofstream *allParamsDown;
  char fileName[200];

  for (unsigned int i = 0; i < Nmasspoints; i++) {
    masspoints[i] = (double)atof(masspointsS[i].c_str());
  }
  
  sprintf(fileName,"text/allParams_%s_%dTeV.txt",nameSample[whichtype].c_str(),LHCsqrts);
  allParams = new ofstream(fileName);
  sprintf(fileName,"text/allParamsUp_oneSyst_%s_%dTeV.txt",nameSample[whichtype].c_str(),LHCsqrts);
  allParamsUp = new ofstream(fileName);
  sprintf(fileName,"text/allParamsDown_oneSyst_%s_%dTeV.txt",nameSample[whichtype].c_str(),LHCsqrts);
  allParamsDown = new ofstream(fileName);

  char thePar[10];
  char equalS[1];
  char dashS[1];
  char pmS[3];
  char theLimit1[10];
  char theLimit2[10];
  float fitted, error;

  // Inizia

  sprintf(fileName,"text/paramsPTOverMCJLST_%s125_%dTeV_Default.txt",nameSample[whichtype].c_str(),LHCsqrts);
  ifstream theFile(fileName);
  
  // first store defaults for 125 (in a random array)
  while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
    
    if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) bbup[0] = fitted;
    if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) Tup[0] = fitted;
    if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) bbdueup[0] = fitted; 
    if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) fexpup[0] = fitted;
    if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) nup[0] = fitted;
    if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) mup[0] = fitted;
    if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) ndue[0] = fitted; 
  }
  
  theFile.close();
  

  for (unsigned int i = 0; i < Nmasspoints; i++) {
  
    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_Default.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts);
    ifstream theFile(fileName);
	
    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
      cout << thePar << " " << fitted << " " << error << endl;
      // limit the error to meaningful values
      // if (fabs(error/fitted) > 0.3) error = 0.3*fitted;
      // if (fabs(error/fitted) < 0.01) error = 0.01*fitted;
      if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) { bb[i] = fitted; bberr[i] = error; bberrup[i] = error*error; bberrdown[i] = error*error; }
      if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) { T[i] = fitted; Terr[i] = error; Terrup[i] = error*error; Terrdown[i] = error*error; }
      if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) { bbdue[i] = fitted; bbdueerr[i] = error; bbdueerrup[i] = error*error; bbdueerrdown[i] = error*error; }
      if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) { fexp[i] = fitted; fexperr[i] = error; fexperrup[i] = error*error; fexperrdown[i] = error*error; }
      if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) { n[i] = fitted; nerr[i] = error; nerrup[i] = error*error; nerrdown[i] = error*error; }
      if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) { m[i] = fitted; merr[i] = error; merrup[i] = error*error; mup[i] = error*error; merrdown[i] = error*error; }
      if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) { ndue[i] = fitted; ndueerr[i] = error; ndueerrup[i] = error*error; ndueerrdown[i] = error*error; }
    }

    theFile.close();

    for (int ss = 0; ss < 5; ss++) {
      if (systSources[whichtype][ss] != "") {

	if (isProfiledSource[whichtype][ss] || i == 2) {  // i.e. m = 125
	  
	  if (systSources[whichtype][ss] == "Mela") 
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_Mela06-10.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts);	    
	  else
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_%s.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts,systSources[whichtype][ss].Data());

	  ifstream theFile(fileName);
	  
          // add errors to the parameters from syst (up and down except MELA)
	  while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	    if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) { bberrup[i] += pow(fitted-bb[i],2); if (systSources[whichtype][ss] != "Mela") bberrdown[i] += pow(fitted-bb[i],2); }
	    if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) { Terrup[i] += pow(fitted-T[i],2); if (systSources[whichtype][ss] != "Mela") Terrdown[i] += pow(fitted-T[i],2);}
	    if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) { bbdueerrup[i] += pow(fitted-bbdue[i],2); if (systSources[whichtype][ss] != "Mela") bbdueerrdown[i] += pow(fitted-bbdue[i],2); }
	    if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) { fexperrup[i] += pow(fitted-fexp[i],2); if (systSources[whichtype][ss] != "Mela") fexperrdown[i] += pow(fitted-fexp[i],2); }
	    if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) { nerrup[i] += pow(fitted-n[i],2); if (systSources[whichtype][ss] != "Mela") nerrdown[i] += pow(fitted-n[i],2); }
	    if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) { merrup[i] += pow(fitted-m[i],2); if (systSources[whichtype][ss] != "Mela") merrdown[i] += pow(fitted-m[i],2); }
	    if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) { ndueerrup[i] += pow(fitted-ndue[i],2); if (systSources[whichtype][ss] != "Mela") ndueerrdown[i] += pow(fitted-ndue[i],2); }
	  }

	  theFile.close();
	  
	  if (systSources[whichtype][ss] == "Mela") {
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_Mela00-13.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts);
	    ifstream theFile(fileName);
	  
	    // add errors to the parameters from syst (down for MELA only)
	    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	      if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) bberrdown[i] += pow(fitted-bb[i],2); 
	      if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) Terrdown[i] += pow(fitted-T[i],2);
	      if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) bbdueerrdown[i] += pow(fitted-bbdue[i],2); 
	      if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) fexperrdown[i] += pow(fitted-fexp[i],2); 
	      if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) nerrdown[i] += pow(fitted-n[i],2); 
	      if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) merrdown[i] += pow(fitted-m[i],2); 
	      if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) ndueerrdown[i] += pow(fitted-ndue[i],2); 
	    }

	  theFile.close();

	  }

	} else {   // i.e. m != 125 or systematics which are NOT profiled
	  	  
	  if (systSources[whichtype][ss] == "Mela") 
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s125_%dTeV_Mela06-10.txt",nameSample[whichtype].c_str(),LHCsqrts);	    
	  else
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s125_%dTeV_%s.txt",nameSample[whichtype].c_str(),LHCsqrts,systSources[whichtype][ss].Data());
	  
	  ifstream theFile(fileName);

	  // add errors to the parameters from relative syst variation observed at 125 (up and down except MELA)
	  while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	    if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) { bberrup[i] += pow((fitted-bbup[0])*bb[i]/bbup[0],2); if (systSources[whichtype][ss] != "Mela") bberrdown[i] += pow((fitted-bbup[0])*bb[i]/bbup[0],2); }
	    if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) { Terrup[i] += pow((fitted-Tup[0])*T[i]/Tup[0],2); if (systSources[whichtype][ss] != "Mela") Terrdown[i] += pow((fitted-Tup[0])*T[i]/Tup[0],2);}
	    if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) { bbdueerrup[i] += pow((fitted-bbdueup[0])*bbdue[i]/bbdueup[0],2); if (systSources[whichtype][ss] != "Mela") bbdueerrdown[i] += pow((fitted-bbdueup[0])*bbdue[i]/bbdueup[0],2); }
	    if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) { fexperrup[i] += pow((fitted-fexpup[0])*fexp[i]/fexpup[0],2); if (systSources[whichtype][ss] != "Mela") fexperrdown[i] += pow((fitted-fexpup[0])*fexp[i]/fexpup[0],2); }
	    if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) { nerrup[i] += pow((fitted-nup[0])*n[i]/nup[0],2); if (systSources[whichtype][ss] != "Mela") nerrdown[i] += pow((fitted-nup[0])*n[i]/nup[0],2); }
	    if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) { merrup[i] += pow((fitted-mup[0])*m[i]/mup[0],2); if (systSources[whichtype][ss] != "Mela") merrdown[i] += pow((fitted-mup[0])*m[i]/mup[0],2); }
	    if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) { ndueerrup[i] += pow((fitted-ndueup[0])*ndue[i]/ndueup[0],2); if (systSources[whichtype][ss] != "Mela") ndueerrdown[i] += pow((fitted-ndueup[0])*ndue[i]/ndueup[0],2); }
	  }

	  theFile.close();
	  
	  if (systSources[whichtype][ss] == "Mela") {
	    sprintf(fileName,"text/paramsPTOverMCJLST_%s%d_%dTeV_Mela00-13.txt",nameSample[whichtype].c_str(),int(masspoints[i]),LHCsqrts);
	    ifstream theFile(fileName);

	    // add errors to the parameters from relative syst variation observed at 125 (down, only MELA)
	    while (theFile >> thePar >> equalS >> fitted >> pmS >> error >> theLimit1 >> dashS >> theLimit2) {
	      if (!strcmp(thePar,"bbup") || !strcmp(thePar,"bb") ) bberrdown[i] += pow((fitted-bbup[0])*bb[i]/bbup[0],2); 
	      if (!strcmp(thePar,"Tup") || !strcmp(thePar,"T") ) Terrdown[i] += pow((fitted-Tup[0])*T[i]/Tup[0],2);
	      if (!strcmp(thePar,"bb2up") || !strcmp(thePar,"bb2") ) bbdueerrdown[i] += pow((fitted-bbdueup[0])*bbdue[i]/bbdueup[0],2); 
	      if (!strcmp(thePar,"fexpup") || !strcmp(thePar,"fexp") ) fexperrdown[i] += pow((fitted-fexpup[0])*fexp[i]/fexpup[0],2); 
	      if (!strcmp(thePar,"nup") || !strcmp(thePar,"n") ) nerrdown[i] += pow((fitted-nup[0])*n[i]/nup[0],2); 
	      if (!strcmp(thePar,"mup") || !strcmp(thePar,"m") ) merrdown[i] += pow((fitted-mup[0])*m[i]/mup[0],2); 
	      if (!strcmp(thePar,"n2up") || !strcmp(thePar,"n2") ) ndueerrdown[i] += pow((fitted-ndueup[0])*ndue[i]/ndueup[0],2); 
	    }

	  theFile.close();

	  }
	}
      }
    }

    bbup[i] = bb[i] + sqrt(bberrup[i]);
    bbdown[i] = bb[i] - sqrt(bberrdown[i]);

    bbdueup[i] = bbdue[i] + sqrt(bbdueerrup[i]);
    bbduedown[i] = bbdue[i] - sqrt(bbdueerrdown[i]);

    nup[i] = n[i] + sqrt(nerrup[i]);
    ndown[i] = n[i] - sqrt(nerrdown[i]);

    mup[i] = m[i] + sqrt(merrup[i]);
    mdown[i] = m[i] - sqrt(merrdown[i]);

    fexpup[i] = fexp[i] + sqrt(fexperrup[i]);
    fexpdown[i] = fexp[i] - sqrt(fexperrdown[i]);

    Tup[i] = T[i] + sqrt(Terrup[i]);
    Tdown[i] = T[i] - sqrt(Terrdown[i]);

    ndueup[i] = ndue[i] + sqrt(ndueerrup[i]);
    nduedown[i] = ndue[i] - sqrt(ndueerrdown[i]);

  }

  TF1 *mypol2 = new TF1("mypol2","pol2");
  TF1 *mypol3 = new TF1("mypol3","pol3");
  TF1 *myZZfunc = new TF1("myZZfunc",zzFunc,100.,1600.,5);
  // TF1 *mypol4 = new TF1("mypol4","pol4");
 
  TCanvas allParCan("allParCan","allParCan",10,10,600,1000);
  allParCan.Divide(2,4);

  allParCan.cd(1);
  fitAndWrite(masspoints,bb,masspointserr,bberr,"bb",myZZfunc,allParams);  

  allParCan.cd(2);
  fitAndWrite(masspoints,T,masspointserr,Terr,"T",myZZfunc,allParams);

  allParCan.cd(3);
  fitAndWrite(masspoints,bbdue,masspointserr,bbdueerr,"bbdue",myZZfunc,allParams);

  allParCan.cd(4);
  fitAndWrite(masspoints,fexp,masspointserr,fexperr,"fexp",myZZfunc,allParams);

  allParCan.cd(5);
  fitAndWrite(masspoints,n,masspointserr,nerr,"n",myZZfunc,allParams);

  allParCan.cd(6);
  fitAndWrite(masspoints,m,masspointserr,merr,"m",myZZfunc,allParams);

  allParCan.cd(7);
  fitAndWrite(masspoints,ndue,masspointserr,ndueerr,"ndue",myZZfunc,allParams);

  sprintf(fileName,"figs/fitParCJLST_%s_%dTeV.pdf",nameSample[whichtype].c_str(),LHCsqrts);
  allParCan.SaveAs(fileName);
  allParams->close();

  allParCan.cd(1);
  fitAndWrite(masspoints,bbup,masspointserr,bberr,"bb",myZZfunc,allParamsUp);  

  allParCan.cd(2);
  fitAndWrite(masspoints,Tup,masspointserr,Terr,"T",myZZfunc,allParamsUp);

  allParCan.cd(3);
  fitAndWrite(masspoints,bbdueup,masspointserr,bbdueerr,"bbdue",myZZfunc,allParamsUp);

  allParCan.cd(4);
  fitAndWrite(masspoints,fexpup,masspointserr,fexperr,"fexp",myZZfunc,allParamsUp);

  allParCan.cd(5);
  fitAndWrite(masspoints,nup,masspointserr,nerr,"n",myZZfunc,allParamsUp);

  allParCan.cd(6);
  fitAndWrite(masspoints,mup,masspointserr,merr,"m",myZZfunc,allParamsUp);

  allParCan.cd(7);
  fitAndWrite(masspoints,ndueup,masspointserr,ndueerr,"ndue",myZZfunc,allParamsUp);

  sprintf(fileName,"figs/fitParCJLST_%s_%dTeV_Up.pdf",nameSample[whichtype].c_str(),LHCsqrts);
  allParCan.SaveAs(fileName);
  allParamsUp->close();

  allParCan.cd(1);
  fitAndWrite(masspoints,bbdown,masspointserr,bberr,"bb",myZZfunc,allParamsDown);  

  allParCan.cd(2);
  fitAndWrite(masspoints,Tdown,masspointserr,Terr,"T",myZZfunc,allParamsDown);

  allParCan.cd(3);
  fitAndWrite(masspoints,bbduedown,masspointserr,bbdueerr,"bbdue",myZZfunc,allParamsDown);

  allParCan.cd(4);
  fitAndWrite(masspoints,fexpdown,masspointserr,fexperr,"fexp",myZZfunc,allParamsDown);

  allParCan.cd(5);
  fitAndWrite(masspoints,ndown,masspointserr,nerr,"n",myZZfunc,allParamsDown);

  allParCan.cd(6);
  fitAndWrite(masspoints,mdown,masspointserr,merr,"m",myZZfunc,allParamsDown);

  allParCan.cd(7);
  fitAndWrite(masspoints,nduedown,masspointserr,ndueerr,"ndue",myZZfunc,allParamsDown);

  sprintf(fileName,"figs/fitParCJLST_%s_%dTeV_Down.pdf",nameSample[whichtype].c_str(),LHCsqrts);
  allParCan.SaveAs(fileName);
  allParamsDown->close();

  return;
}
