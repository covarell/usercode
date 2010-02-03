#include <iostream>
#include <string.h>
#include <iomanip>
#include<fstream>
#include <math.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"


void effPlotsData(string runNumber){

  string MCorRECO = "RECO";  
  string f1 = "DQM_V0001_HLT_R" ;
  int howManyZeros = 9 - strlen(runNumber.c_str());
  for (int i=0;i<howManyZeros;i++) f1 = f1 + "0";
  f1 = f1 + runNumber + ".root";
  cout << "Input file: " << f1 << endl;

  TFile* file1 = TFile::Open(f1.c_str());
  // TFile* file2 = TFile::Open(f2);

  TH1F* hf1=0, *hf2=0;

  TH1F* num1=0, *denom1=0, *num2=0, *denom2=0; 

  char* name1=0, *name2=0;

  const int npaths = 2;
  string paths[npaths]={"HLT_L1SingleEG2_DQMGammaJet","HLT_Photon10_L1R_DQMGammaJet"};
  string Histo_paths[npaths];
  for(int i=0;i<npaths;i++){
    Histo_paths[i]="DQMData/Run " + runNumber + "/HLT/Run summary/HLTEgammaValidation/"+paths[i]+"/";
    cout << "Path " << i << ": " << Histo_paths[i] << endl;
  }

  string outname;
  outname = "effPlots" + runNumber + ".root"; 
  TFile fout(outname.c_str(),"RECREATE");
  TObjArray histarr(0);
  cout << "Writing to file: " << outname << endl;
  //TCanvas* c1= new TCanvas("c1","c1");
  //c1->SetFillColor(0);
  // c1->Divide(2,2);

  for(int path=0; path< npaths ; path++){
    
    string firstname = Histo_paths[path] + "total_eff_" + MCorRECO + "_matched";
    hf1 = (TH1F*)file1->Get(firstname.c_str());
    
    name2 = hf1->GetXaxis()->GetBinLabel(1);
    std::cout << name2 << std::endl;   
    
    string numname = Histo_paths[path] + name2 + "et_" + MCorRECO + "_matched";
    string denomname;
    if (MCorRECO == "RECO") denomname = Histo_paths[path] + "reco_et";
    else denomname = Histo_paths[path] + "gen_et";
    string numnameeta = Histo_paths[path] + name2 + "eta_" + MCorRECO + "_matched";
    string denomnameeta;
    if (MCorRECO == "RECO") denomnameeta = Histo_paths[path] + "reco_eta";
    else denomnameeta = Histo_paths[path] + "gen_eta";

    std::cout << denomname << std::endl;  

    num1   =  new TH1F( *(TH1F*)file1->Get(numname.c_str()) );
    denom1 =  new TH1F( *(TH1F*)file1->Get(denomname.c_str()) );
    num2   =  new TH1F( *(TH1F*)file1->Get(numnameeta.c_str()) );
    denom2 =  new TH1F( *(TH1F*)file1->Get(denomnameeta.c_str()) );
    // num2   =  new TH1F( *(TH1F*)file2->Get(numname.c_str()) );
    // std::cout << name1 << std::endl;  
    
    num1->Sumw2();
    denom1->Sumw2();
    num2->Sumw2();
    denom2->Sumw2();
    
    num1->Divide(num1,denom1,1.,1.,"b");
    num2->Divide(num2,denom2,1.,1.,"b");
    
    string title = paths[path] + " : " + name2 + "; et ; Efficiency";
    string thename = paths[path] + name2 + "eteff";
    num1->SetName(thename.c_str());
    num1->SetTitle(title.c_str());
    num1->SetMarkerColor(2);
    num1->SetMarkerStyle(20);
    title = paths[path] + " : " + name2 + "; eta ; Efficiency";
    thename = paths[path] + name2 + "etaeff";
    num2->SetName(thename.c_str());
    num2->SetTitle(title.c_str());
    num2->SetMarkerColor(2);
    num2->SetMarkerStyle(20);
    
    histarr.Add(num1);
    histarr.Add(num2);

    for(int filter=1; filter < hf1->GetNbinsX()-2 ; filter++){
      name1 = hf1->GetXaxis()->GetBinLabel(filter);
      name2 = hf1->GetXaxis()->GetBinLabel(filter+1);
      std::cout << name1 << std::endl;   
      
      string numname = Histo_paths[path] + name2 + "et_" + MCorRECO + "_matched";
      string denomname = Histo_paths[path] + name1 + "et_" + MCorRECO + "_matched";
      string numnameeta = Histo_paths[path] + name2 + "eta_" + MCorRECO + "_matched";
      string denomnameeta = Histo_paths[path] + name1 + "eta_" + MCorRECO + "_matched";
      
      std::cout << denomname << std::endl;  

      num1   =  new TH1F( *(TH1F*)file1->Get(numname.c_str()) );
      denom1 =  new TH1F( *(TH1F*)file1->Get(denomname.c_str()) );
      num2   =  new TH1F( *(TH1F*)file1->Get(numnameeta.c_str()) );
      denom2 =  new TH1F( *(TH1F*)file1->Get(denomnameeta.c_str()) );
      // num2   =  new TH1F( *(TH1F*)file2->Get(numname.c_str()) );
      // std::cout << name1 << std::endl;  
      
      num1->Sumw2();
      denom1->Sumw2();
      num2->Sumw2();
      denom2->Sumw2();

      num1->Divide(num1,denom1,1.,1.,"b");
      num2->Divide(num2,denom2,1.,1.,"b");
   
      string title = paths[path] + " : " + name2 + "; et ; Efficiency";
      string thename = paths[path] + name2 + "eteff";
      num1->SetName(thename.c_str());
      num1->SetTitle(title.c_str());
      num1->SetMarkerColor(2);
      num1->SetMarkerStyle(20);
      title = paths[path] + " : " + name2 + "; eta ; Efficiency";
      thename = paths[path] + name2 + "etaeff";
      num2->SetName(thename.c_str());
      num2->SetTitle(title.c_str());
      num2->SetMarkerColor(2);
      num2->SetMarkerStyle(20);
 
      histarr.Add(num1);
      histarr.Add(num2);
      /* TLegend *legend = new TLegend(0.4,0.2,0.55,0.4,"");
      TLegendEntry* l1 = legend->AddEntry(num1,label1,"l");
      l1->SetTextSize(0.1);
      l1 = legend->AddEntry(num2,label2,"l");
      l1->SetTextSize(0.1);*/
    }
  }
  fout.cd();
  histarr.Write();
  fout.Close();
}


