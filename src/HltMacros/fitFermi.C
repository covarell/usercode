#include <iostream>
#include <string.h>
#include <iomanip>
#include<fstream>
#include <math.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"


void fitFermi(char* f1, string histoname){
  TFile* file1 = TFile::Open(f1);
  
  TH1F* hf1 = (TH1F*)file1->Get(histoname.c_str());
  hf1->SetMarkerColor(4);
  for (int i=1;i<hf1->GetNbinsX();i++) {
    if (hf1->GetBinContent(i) == 1) hf1->SetBinError(i,0.5);
  }
  TF1* fermif = new TF1("fermif","[0]/(1+exp(([1]-x)/[2]))",0.,20.);
  fermif->SetParameter(0,1.);
  fermif->SetParameter(1,4.0);
  fermif->SetParameter(2,0.5);

  hf1->Fit("fermif","","e",1.,18.);   
  //fermif->Draw();
  
}


