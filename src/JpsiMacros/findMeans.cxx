// ROOT includes
#include <TROOT.h>
#include <TFile.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>

void getrange(string &varRange, float *varmin, float *varmax)
{
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
   cout << varRange.c_str() << ": range not valid!" << endl;
    assert(0);
  }

 return;
}

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");

  char *filename;
  string prange;
  string etarange;
  int theTrigger = -1;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];
    
    switch(pchar[0]){
      
    case '-':{
      
      switch(pchar[1]){
	
      case 'f':{
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;
      }
     
      case 'p':{
	prange = argv[i+1];
	cout << "Range for pT is " << prange << " GeV/c" << endl;
        break;
      }
       
      case 'e':{
        etarange = argv[i+1];
        cout << "Range for |eta| is " << etarange << endl;
        break;
      }

      case 't':{
        theTrigger = atoi(argv[i+1]);
        cout << "Using trigger bit n. " << theTrigger << endl;
        break; 

      }
      }
    }
    }
  }
  
  TFile fIn(filename);
  fIn.cd();
  
  RooDataSet *data = (RooDataSet*)fIn.Get("data");
    
  float pmin, pmax; 
  float etamin, etamax;
  
  getrange(prange,&pmin,&pmax);
  getrange(etarange,&etamin,&etamax);

  RooDataSet *reddata;
 
  char reducestr[200];
 
  char oFile[200];
  sprintf(oFile,"results/meanPt/results_pT%s_eta%s.txt",prange.c_str(),etarange.c_str());
  ofstream outputFile(oFile);
  
  // for (int j=0; j < 2; j++) {

  sprintf(reducestr,"JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f", pmax,pmin,etamax,etamin);
  
  // for selecting triggers and only opposite sign pairs if needed 
  const RooArgSet* thisRow = data->get(0);  
  RooCategory* theSign = (RooCategory*)thisRow->find("JpsiSign");
  if (theSign) sprintf(reducestr,"%s && JpsiSign == JpsiSign::OS",reducestr);
  RooCategory* thetrigger = (RooCategory*)thisRow->find("triggerMu");
  if (thetrigger) {
    if (theTrigger == 0) sprintf(reducestr,"%s && triggerDMu > 0",reducestr); 
    else if (theTrigger == 1) sprintf(reducestr,"%s && triggerMuPre > 0",reducestr); 
    else if (theTrigger == 2) sprintf(reducestr,"%s && triggerMu > 0",reducestr);
    else if (theTrigger == 3) sprintf(reducestr,"%s && triggerOniaTrack > 0",reducestr);
    else if (theTrigger == 4) sprintf(reducestr,"%s && triggerOniaL1Mu > 0",reducestr);
  }
  //
  
  reddata = (RooDataSet*)data->reduce(reducestr);
  // reddata->setWeightVar("MCweight");
  
  float meanP = 0.;
  float meanP2 = 0.;
  
  // const RooArgSet* thisRow;   
  
  for (Int_t iSamp = 0; iSamp < reddata->numEntries(); iSamp++)
    {
      thisRow = reddata->get(iSamp);
      
      RooRealVar* myPt = (RooRealVar*)thisRow->find("JpsiPt");
      float thePt = myPt->getVal();
      float theWeight = reddata->weight();
      meanP += thePt*theWeight;
      meanP2 += thePt*thePt*theWeight;
      
    }
  
  // reddata->setWeightVar("MCweight");
  
  if (reddata->sumEntries() > 0 ) {
    meanP = meanP / reddata->sumEntries();
    meanP2 = meanP2 / reddata->sumEntries();
  }
  
  float rms = sqrt(meanP2 - pow(meanP,2));
  outputFile << "ME " << meanP << " ME2 " << meanP2 << " RMS " << rms << endl;

  return 0;
  
}

