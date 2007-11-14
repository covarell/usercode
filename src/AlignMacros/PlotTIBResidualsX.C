void PlotTIBResidualsX() {

  // MY SETTINGS
  TString fileName("/afs/cern.ch/user/c/covarell/scratch0/joboutput/Pass3TIBmod-x-outrej-linAPE/CMSSW_1_3_6/main/CSA06AlignmentEvents.root");
  int theIterations[3] = {1,3,10};
  float plotRange = 0.4; 

  TFile f(fileName.Data());
  char theTreeName[20];
  char theHistoName[20];

  h1_1 = new TH1F("h1_1","Residuals",80,-plotRange,plotRange);
  h1_1->SetLineWidth(2);
  h2_1 = new TH1F("h2_1","Residuals",80,-plotRange,plotRange);
  h2_1->SetLineColor(kRed); 
  h2_1->SetLineWidth(2);
  h3_1 = new TH1F("h3_1","Residuals",80,-plotRange,plotRange);
  h3_1->SetLineColor(kBlue); 
  h3_1->SetLineWidth(2);
  h1_2 = new TH1F("h1_2","Residuals",80,-plotRange,plotRange);
  h1_2->SetLineWidth(2);
  h2_2 = new TH1F("h2_2","Residuals",80,-plotRange,plotRange);
  h2_2->SetLineColor(kRed); 
  h2_2->SetLineWidth(2);
  h3_2 = new TH1F("h3_2","Residuals",80,-plotRange,plotRange);
  h3_2->SetLineColor(kBlue); 
  h3_2->SetLineWidth(2);
  h1_3 = new TH1F("h1_3","Residuals",80,-plotRange,plotRange);
  h1_3->SetLineWidth(2);
  h2_3 = new TH1F("h2_3","Residuals",80,-plotRange,plotRange);
  h2_3->SetLineColor(kRed); 
  h2_3->SetLineWidth(2);
  h3_3 = new TH1F("h3_3","Residuals",80,-plotRange,plotRange);
  h3_3->SetLineColor(kBlue); 
  h3_3->SetLineWidth(2);
  h1_4 = new TH1F("h1_4","Residuals",80,-plotRange,plotRange);
  h1_4->SetLineWidth(2);
  h2_4 = new TH1F("h2_4","Residuals",80,-plotRange,plotRange);
  h2_4->SetLineColor(kRed); 
  h2_4->SetLineWidth(2);
  h3_4 = new TH1F("h3_4","Residuals",80,-plotRange,plotRange);
  h3_4->SetLineColor(kBlue); 
  h3_4->SetLineWidth(2);

  sprintf(theTreeName,"T1_%d",theIterations[0]);
  TTree * tree1 = (TTree *) gROOT->FindObject(theTreeName);
  sprintf(theTreeName,"T1_%d",theIterations[1]);
  TTree * tree2 = (TTree *) gROOT->FindObject(theTreeName);
  sprintf(theTreeName,"T1_%d",theIterations[2]);
  TTree * tree3 = (TTree *) gROOT->FindObject(theTreeName);

  TCanvas c1("c1","Residuals",10,10,1100,800);
  c1.Divide(2,2);

  char theLayer[100];
  
  for (int i = 1; i <= 4; i++) {

    // ** Do not apply cuts
    // sprintf(theLayer,"hLayer == %d && hType == 3",i);
    // ** Apply cuts
    sprintf(theLayer,"hLayer == %d && hType == 3 && isOnAli == 1",i);
    sprintf(theHistoName,"xres >> h3_%d",i);
    tree3->Draw(theHistoName,theLayer,"goff");   
    sprintf(theHistoName,"xres >> h2_%d",i);
    tree2->Draw(theHistoName,theLayer,"goff");   
    sprintf(theHistoName,"xres >> h1_%d",i);
    tree1->Draw(theHistoName,theLayer,"goff");   
  }

  // Normalize
  /* h1_1->Scale(1./h1_1->Integral());
  h2_1->Scale(1./h2_1->Integral());
  h3_1->Scale(1./h3_1->Integral());
  h1_2->Scale(1./h1_2->Integral());
  h2_2->Scale(1./h2_2->Integral());
  h3_2->Scale(1./h3_2->Integral());
  h1_3->Scale(1./h1_3->Integral());
  h2_3->Scale(1./h2_3->Integral());
  h3_3->Scale(1./h3_3->Integral());
  h1_4->Scale(1./h1_4->Integral());
  h2_4->Scale(1./h2_4->Integral());
  h3_4->Scale(1./h3_4->Integral()); */

  c1.cd(1);
  h3_1->Draw();
  h2_1->Draw("SAME");
  h1_1->Draw("SAME");
  c1.cd(2);
  h3_2->Draw();
  h2_2->Draw("SAME");
  h1_2->Draw("SAME");
  c1.cd(3);
  h3_3->Draw();
  h2_3->Draw("SAME");
  h1_3->Draw("SAME");
  c1.cd(4);
  h3_4->Draw();
  h2_4->Draw("SAME");
  h1_4->Draw("SAME");

  c1.Print("TIBresidualsX.eps");
}
