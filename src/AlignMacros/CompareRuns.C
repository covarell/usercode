#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <fstream.h>

#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TPostScript.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TPaveLabel.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom.h>

Int_t CompareRuns();

#endif

#include "IOAlignableDataTree.h"
#include "IOPrivateTree.h"
#include "IOCombined.h"

// ----------------------------------------------------------------------------

Int_t CompareRuns() {

   // make nice plots ---------------------------------------------------------

  // general style for all plots
   gStyle->SetOptStat(1111111);    // contents of statistics box
   gStyle->SetHistLineWidth(2); // histo line width
   //gROOT->SetStyle("Plain");    // no grey background 

   // my style for 2x3 plots
   TStyle* style23  = new TStyle("style23","My style for 2x3"); 
   style23->SetCanvasBorderMode(0);
   style23->SetCanvasColor(0);

   style23->SetPadBorderMode(0);
   style23->SetPadBorderSize(0);
   style23->SetPadColor(0);

   // for 1x1 or 2x2
   // axes labels
   //style23->SetLabelSize(0.04,"X"); 
   //style23->SetLabelSize(0.06,"Y");
   // axes title
   //style23->SetTitleSize(0.045,"X");
   //style23->SetTitleSize(0.045,"Y");

   // for 2x3
   // axes labels
   style23->SetLabelSize(0.075,"X"); 
   style23->SetLabelSize(0.075,"Y");
   // axes title
   style23->SetTitleSize(0.085,"YZ");
   style23->SetTitleSize(0.085,"X");
   //style23->SetTitleSize(0.065,"Y");

   style23->SetTitleXOffset(0.1);
   style23->SetTitleYOffset(0.8);

  style23->SetPadTopMargin(0.05);
  style23->SetPadBottomMargin(0.2); // 0.13
  style23->SetPadLeftMargin(0.16);
  style23->SetPadRightMargin(0.02);

   // title
   style23->SetTitleFontSize(0.07);
   style23->SetTitleOffset(0.);
   style23->SetTitleStyle(1001);
   style23->SetTitleBorderSize(1);  
   style23->SetTitleColor(1);
   style23->SetTitleFillColor(10);  

   style23->SetTextSize(0.075);

   style23->SetStatBorderSize(1);  
   style23->SetStatColor(0);  
   style23->SetOptStat(0);    // contents of statistics box
   style23->SetOptFit(11110); // contents of fit box
   style23->SetStatFontSize(0.075);  
   style23->SetStatX(0.95);
   style23->SetStatY(0.92);

   style23->SetHistLineWidth(2); // histo line width

   style23->SetPalette(1);

   // set style
   gROOT->SetStyle("style23");

   IOCombined *ioc1 = new IOCombined();   
   IOCombined *ioc2 = new IOCombined(); 
   IOCombined *ioc3 = new IOCombined();   
   IOCombined *ioc4 = new IOCombined(); 

   TString psname = "psfiles/comparison-confA-confC.ps";
   TH1F *halidiffxTOB = new TH1F("halidiffxTOB","; #Deltax (runs1 - runs2) [#mum]", 10, -300.,300.);
   TH1F *halidiffyTOB = new TH1F("halidiffyTOB","; #Deltay (runs1 - runs2) [#mum]", 10, -900.,1200.); 
   TH1F *halidiffxTIB = new TH1F("halidiffxTIB","; #Deltax (runs1 - runs2) [#mum]", 10, -300.,300.);
   TH1F *halidiffyTIB = new TH1F("halidiffyTIB","; #Deltay (runs1 - runs2) [#mum]", 10, -900.,1200.);  
   // TH1F *halidiffz = new TH1F("halidiffz","; #Deltaz (runs1 - runs2) [#mum]", 10, -500.,500.); 
   // TH1F *halidiffb = new TH1F("halidiffb","; #Delta#beta (runs1 - runs2) [rad]", 10, -0.0005,0.0005); 
   // TH1F *halidiffc = new TH1F("halidiffc","; #Delta#gamma (all - 6217) [rad]", 20, -0.04,0.04);
   
   // read in and process input trees - TOB --------------------------
   // FIRST SET
   TString dir3 = "/afs/cern.ch/user/c/covarell/scratch0/joboutput/Pass3TOBmod-x-outrej-linAPE/CMSSW_1_3_6/main/";
   TString dir4 = "/afs/cern.ch/user/c/covarell/scratch0/joboutput/Pass3CTOBmod-xy-outrej-linAPE/CMSSW_1_3_6/main/";
   TString nametree3 = "AlignablesAbsPos_10"; 
   TString nametree4 = "AlignablesAbsPos_10";
   
   TFile* filetrue3 = new TFile(dir3+"IOTruePositions.root");
   TTree* treetrue3 = (TTree*)filetrue3->Get("AlignablesOrgPos_1");
   if (treetrue3 !=0) {
     printf("filling true positions 1 ...\n");
     IOAlignableDataTree* alitrue3 = new IOAlignableDataTree(treetrue3);
       alitrue3->fill(ioc3,1,0); // 1 means true positions 
   }
   else { cout << "Tree not found!\n"; return 0; } 

   TFile* file3 = new TFile(dir3+"IOAlignedPositions.root");
   TTree* treeAligned3 = (TTree*)file3->Get(nametree3);
   if (treeAligned3 !=0) {
     printf("filling aligned positions 1 ...\n");
     IOAlignableDataTree* alitree3 = new IOAlignableDataTree(treeAligned3);
       alitree3->fill(ioc3,3,1); // 3 means aligned positions, 1 is last iter.
   }
   else { cout << "Tree not found!\n"; return 0; }
  
   // SECOND SET

   TFile* filetrue4 = new TFile(dir4+"IOTruePositions.root");
   TTree* treetrue4 = (TTree*)filetrue4->Get("AlignablesOrgPos_1");
   if (treetrue4 !=0) {
     printf("filling true positions 2 ...\n");
     IOAlignableDataTree* alitrue4 = new IOAlignableDataTree(treetrue4);
       alitrue4->fill(ioc4,1,0); // 1 means true positions 
   }
   else { cout << "Tree not found!\n"; return 0; } 

   TFile* file4 = new TFile(dir4+"IOAlignedPositions.root");
   TTree* treeAligned4 = (TTree*)file4->Get(nametree4);
   if (treeAligned4 !=0) {
     printf("filling aligned positions 2 ...\n");
     IOAlignableDataTree* alitree4 = new IOAlignableDataTree(treeAligned4);
       alitree4->fill(ioc4,3,1); // 3 means aligned positions, 1 is last iter.
   }
   else { cout << "Tree not found!\n"; return 0; }

   printf("filling private tree ...\n");
   TFile* filePriv2 = new TFile(dir4+"CSA06AlignmentAlignables.root");
   TTree* treePriv2 = (TTree*)filePriv2->Get("T2");
   if (treePriv2==0) {
     cout <<"Private tree not found!\n";
     return 0;
   }
   IOPrivateTree* privtree2 = new IOPrivateTree(treePriv2);
   privtree2->fill(ioc4);

   for (Int_t i=0; i<ioc3->NAli; i++) {
     for (Int_t j=0; j<ioc4->NAli; j++) {
       if ( (ioc3->AliId[i] == ioc4->AliId[j]) && (ioc3->AliObjId[i] == ioc4->AliObjId[j])) {
         float thedx = 10000*(ioc3->AbsPosAl[i][1][0] - ioc4->AbsPosAl[j][1][0]);
	 float thedy = 10000*(ioc3->AbsPosAl[i][1][1] - ioc4->AbsPosAl[j][1][1]);
	 // float thedz = 10000*(ioc1->AbsPosAl[i][1][2] - ioc2->AbsPosAl[j][1][2]);
         // float thedb = ioc1->AbsPosAl[i][1][4] - ioc2->AbsPosAl[j][1][4];
         cout << j << " " << 
	   ioc4->AliType[j] << " " << 
	   ioc4->AliLayer[j] << " " << 
	   thedx << " " << thedy << " " << endl;
	 if (ioc4->AliType[j] == 5) {
	   if (fabs(thedx) > 0.001)  halidiffxTOB->Fill(thedx);
	   if (fabs(thedy) > 0.001 && ioc4->AliLayer[j] < 3) halidiffyTOB->Fill(thedy);
	 }
         // if (fabs(thedz) > 0.1) halidiffz->Fill(thedz);
         // if (fabs(thedb) > 0.0001) halidiffb->Fill(thedb);
        
       }
     }
   }

  // read in and process input trees - TIB ------------------------------------
   
   // FIRST SET
   TString dir1 = "/afs/cern.ch/user/c/covarell/scratch0/joboutput/Pass3TIBmod-x-outrej-20iter/CMSSW_1_3_6/main/";
   TString dir2 = "/afs/cern.ch/user/c/covarell/scratch0/joboutput/Pass3CTIBmod-xy-outrej-20iter/CMSSW_1_3_6/main/";
   TString nametree1 = "AlignablesAbsPos_20"; 
   TString nametree2 = "AlignablesAbsPos_19";
   
   TFile* filetrue1 = new TFile(dir1+"IOTruePositions.root");
   TTree* treetrue1 = (TTree*)filetrue1->Get("AlignablesOrgPos_1");
   if (treetrue1 !=0) {
     printf("filling true positions 1 ...\n");
     IOAlignableDataTree* alitrue1 = new IOAlignableDataTree(treetrue1);
       alitrue1->fill(ioc1,1,0); // 1 means true positions 
   }
   else { cout << "Tree not found!\n"; return 0; } 

   TFile* file1 = new TFile(dir1+"IOAlignedPositions.root");
   TTree* treeAligned1 = (TTree*)file1->Get(nametree1);
   if (treeAligned1 !=0) {
     printf("filling aligned positions 1 ...\n");
     IOAlignableDataTree* alitree1 = new IOAlignableDataTree(treeAligned1);
       alitree1->fill(ioc1,3,1); // 3 means aligned positions, 1 is last iter.
   }
   else { cout << "Tree not found!\n"; return 0; }
  
   // SECOND SET

   TFile* filetrue2 = new TFile(dir2+"IOTruePositions.root");
   TTree* treetrue2 = (TTree*)filetrue2->Get("AlignablesOrgPos_1");
   if (treetrue2 !=0) {
     printf("filling true positions 2 ...\n");
     IOAlignableDataTree* alitrue2 = new IOAlignableDataTree(treetrue2);
       alitrue2->fill(ioc2,1,0); // 1 means true positions 
   }
   else { cout << "Tree not found!\n"; return 0; } 

   TFile* file2 = new TFile(dir2+"IOAlignedPositions.root");
   TTree* treeAligned2 = (TTree*)file2->Get(nametree2);
   if (treeAligned2 !=0) {
     printf("filling aligned positions 2 ...\n");
     IOAlignableDataTree* alitree2 = new IOAlignableDataTree(treeAligned2);
       alitree2->fill(ioc2,3,1); // 3 means aligned positions, 1 is last iter.
   }
   else { cout << "Tree not found!\n"; return 0; }

   printf("filling private tree ...\n");
   TFile* filePriv = new TFile(dir2+"CSA06AlignmentAlignables.root");
   TTree* treePriv = (TTree*)filePriv->Get("T2");
   if (treePriv==0) {
     cout <<"Private tree not found!\n";
     return 0;
   }
   IOPrivateTree* privtree = new IOPrivateTree(treePriv);
   privtree->fill(ioc2);

   for (Int_t i=0; i<ioc1->NAli; i++) {
     for (Int_t j=0; j<ioc2->NAli; j++) {
       if ( (ioc1->AliId[i] == ioc2->AliId[j]) && (ioc1->AliObjId[i] == ioc2->AliObjId[j])) {
         float thedx = 10000*(ioc1->AbsPosAl[i][1][0] - ioc2->AbsPosAl[j][1][0]);
	 float thedy = 10000*(ioc1->AbsPosAl[i][1][1] - ioc2->AbsPosAl[j][1][1]);
	 // float thedz = 10000*(ioc1->AbsPosAl[i][1][2] - ioc2->AbsPosAl[j][1][2]);
         // float thedb = ioc1->AbsPosAl[i][1][4] - ioc2->AbsPosAl[j][1][4];
         cout << j << " " << 
	   ioc2->AliType[j] << " " << 
	   ioc2->AliLayer[j] << " " << 
	   thedx << " " << thedy << " " << endl;
	 if (ioc2->AliType[j] == 3) {
	   if (fabs(thedx) > 0.001)  halidiffxTIB->Fill(thedx);
	   if (fabs(thedy) > 0.001 && ioc2->AliLayer[j] < 3) halidiffyTIB->Fill(thedy);
	 }
         // if (fabs(thedz) > 0.1) halidiffz->Fill(thedz);
         // if (fabs(thedb) > 0.0001) halidiffb->Fill(thedb);
        
       }
     }
   }

  // open ps output file
  gSystem->Exec("rm -f "+psname);
  TCanvas *can = new TCanvas("can", "Alignment");

  can->Clear();
  can->Divide(2,2);
  can->cd(1); 
  // halidiffx->Draw();
  halidiffxTOB->Fit("gaus","","",-250.,250.);
  can->cd(2); 
  // halidiffy->Draw();
  halidiffyTOB->Fit("gaus","","",-700.,700.);
  can->cd(3); 
  // halidiffx->Draw();
  halidiffxTIB->Fit("gaus","","",-250.,250.);
  can->cd(4); 
  // halidiffy->Draw();
  halidiffyTIB->Fit("gaus","","",-700.,700.);

  can->Print(psname); 
  printf("Finished.\n");
  return 0;
}



