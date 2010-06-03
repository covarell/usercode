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
 
  string type = "GG";
  char oFile[200];
  sprintf(oFile,"results/meanPt/results_pT%s_eta%s.txt",prange.c_str(),etarange.c_str());
  ofstream outputFile(oFile);
  
  for (int j=0; j < 2; j++) {

    if (j == 1) type = "GT";

    sprintf(reducestr,"JpsiPt < %f && JpsiPt > %f && abs(JpsiEta) < %f && abs(JpsiEta) > %f && JpsiType == JpsiType::%s", pmax,pmin,etamax,etamin,type.c_str());

    reddata = (RooDataSet*)data->reduce(reducestr);
    // reddata->setWeightVar("MCweight");
    
    float meanP = 0.;
    float meanP2 = 0.;
  
    const RooArgSet* thisRow;   
  
    for (Int_t iSamp = 0; iSamp < reddata->numEntries(); iSamp++)
      {
	thisRow = reddata->get(iSamp);

	RooRealVar* myPt = (RooRealVar*)thisRow->find("JpsiPt");
        float thePt = myPt->getVal();
        float theWeight = reddata->weight();
	meanP += thePt*theWeight;
	meanP2 += thePt*thePt*theWeight;
        
      }
    
    reddata->setWeightVar("MCweight");

    if (reddata->sumEntries() > 0 ) {
       meanP = meanP / reddata->sumEntries();
       meanP2 = meanP2 / reddata->sumEntries();
    }

    outputFile << type << " " << meanP << " " << meanP2 << " 0.0" << endl;
    
  }

  return 0;
  
}

