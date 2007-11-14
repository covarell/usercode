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

Int_t Alignment();

#endif

void PlotTLHisto(TH2F* histo, Option_t* opt);
void PlotTLHisto(TH1F* histo);

#include "IOAlignableDataTree.h"
#include "IOPrivateTree.h"
#include "IOEventTree.h"
#include "IOAlignmentParameterTree.h"
#include "IOCombined.h"

#include "Analysis.h"

// ----------------------------------------------------------------------------

Int_t Alignment() {

   // make nice plots ---------------------------------------------------------

  // general style for all plots
   gStyle->SetOptStat(1110);    // contents of statistics box
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
   style23->SetOptStat(1110);    // contents of statistics box
   style23->SetStatFontSize(0.075);  
   style23->SetStatX(0.95);
   style23->SetStatY(0.92);


   style23->SetHistLineWidth(2); // histo line width

   style23->SetPalette(1);


   // set style
   gROOT->SetStyle("style23");

   // read in and process input trees -----------------------------------------

   TString path="/afs/cern.ch/user/c/covarell/scratch0/joboutput/";
   TString orca="/CMSSW_1_3_6/";
   orca+="main/";
   
   // DEFAULT VALUES (DO NOT CHANGE) ------------------------------------------
   Int_t MinHit=100;
   Int_t IterFile=50;   // default value, then read from job output
        
   Bool_t PxZoom=kFALSE;
   Float_t ms[2],mr[2];
   //max y for the plot shift 
   mr[0]=0.00149; mr[1]=0.00149 ; ms[0]=350; ; ms[1]=350;

   Int_t itp[4];
   itp[0]=0; itp[1]=5; itp[2]=10; itp[3]=15; // itp[4]=15;
   Int_t ioff=0;
   Int_t nbin=30;

 //  --------------------------------------------------------------------------
 //  MY SETTINGS --------------------------------------------------------------
 //  --------------------------------------------------------------------------

   TString baseJob = "Pass3TIBTOBstr";
   TString floatVar = "xy";    //can contain "xyzabg" in any order
   TString mySetting = "allfree-linAPE";  

   nbin=51;
   PxZoom=false;  
   MinHit=250; 
   ms[0]=500; 
   ms[1]=500; 
   mr[0]=0.0005; 
   mr[1]=0.0005;
   itp[0]=0; itp[1]=1; itp[2]=2; itp[3]=10; 
   // itp[0]=0; itp[1]=5; itp[2]=10; itp[3]=20; 
   
   // Base file name
   TString marker = "-";
   TString job = baseJob+marker+floatVar;
   if (mySetting != "") {
     job += marker;
     job += mySetting;
   } 
 
   // number of iteration (modifiable by hand) 
   TString siterat        =path+job+orca+"IOIteration.root";
   cout << siterat << endl;
   ifstream *numIterat = new ifstream(siterat.Data());
   Int_t aaa = 0;
   if (numIterat) {
     while (!aaa) {
       *numIterat >> IterFile;    aaa++;
     }
   }
 
   // IterFile = 9;
   cout << "Number of iterations found in RootFiles = " << IterFile << endl;

   /////////////----------------------------------------------//////////////////

   //   Open output ps and root files to store all histos 
   char fileNames[50];
   sprintf(fileNames,"rootfiles/%s.root",job.Data()) ;
   // sprintf(fileNames,"%s.root",job.Data()) ;
   TFile * rootFile = new TFile(fileNames, "RECREATE") ;
   sprintf(fileNames,"psfiles/%s.ps",job.Data()) ;
   // sprintf(fileNames,"%s.ps",job.Data()) ; 
   TString psname = fileNames;

// ----------------------------------------------------------------------------

   // TString sfileevt       =path+job+orca+"HIPAlignmentEvents.root";
   TString sfileevt       =path+job+orca+"CSA06AlignmentEvents.root";
   TString sfiletrue      =path+job+orca+"IOTruePositions.root";
   // TString sfiletrue      =path+job+orca+truefile;
   TString sfilemisaligned=path+job+orca+"IOMisalignedPositions.root";
   TString sfilealigned   =path+job+orca+"IOAlignedPositions.root";
   TString sfileparams    =path+job+orca+"IOAlignmentParameters.root";
   // TString sfilepriv      =path+job+orca+"HIPAlignmentAlignables.root";
   TString sfilepriv      =path+job+orca+"CSA06AlignmentAlignables.root";


// end steering - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// create combined object
   IOCombined* ioc=new IOCombined();

// fill true positions
   TFile* fileTrue = new TFile(sfiletrue);
   TTree* treeTrueT1 = (TTree*)fileTrue->Get("AlignablesOrgPos_1");
   IOAlignableDataTree* truetree = new IOAlignableDataTree(treeTrueT1);
   printf("filling true positions ...\n");
   truetree->fill(ioc,1,0); // 1 means true positions
// fill misaligned positions
   TFile* fileMisaligned = new TFile(sfilemisaligned);
   TTree* treeMisalignedT1 = (TTree*)fileMisaligned->Get("AlignablesAbsPos_1");
   IOAlignableDataTree* misalitree = new IOAlignableDataTree(treeMisalignedT1);
   printf("filling misaligned positions ...\n");
   misalitree->fill(ioc,2,0); // 2 means misaligned positions

// fill aligned positions
   TFile* fileAligned = new TFile(sfilealigned);
   //printf("filling aligned positions ...\n");
   for (Int_t i=1+ioff; i<=IterFile; i++) {
     char stree[99];
     sprintf(stree,"AlignablesAbsPos_%d",i);
     TTree* treeAligned = (TTree*)fileAligned->Get(stree);
     if (treeAligned !=0) {
       printf("filling aligned positions %i ...\n",i);
       IOAlignableDataTree* alitree = new IOAlignableDataTree(treeAligned);
       alitree->fill(ioc,3,i-ioff); // 3 means aligned positions, i is iter.
     }
     else { cout << stree <<" not found!\n"; return 0; }
   }

// fill alignment parameters
   TFile* fileParams = new TFile(sfileparams);
   printf("filling alignment parameters ...\n");
   for (Int_t i=1+ioff; i<=IterFile;i++) {
     char stree[99];
     sprintf(stree,"AlignmentParameters_%d",i);
     TTree* treeParams = (TTree*)fileParams->Get(stree);
     if (treeParams !=0) {
       IOAlignmentParameterTree* paratree = 
         new IOAlignmentParameterTree(treeParams);
       paratree->fill(ioc,i-ioff); 
     }
     else { cout <<"not found!\n"; return 0; }
   }

// now fill information from private tree
   printf("filling private tree ...\n");
   TFile* filePriv = new TFile(sfilepriv);
   TTree* treePriv = (TTree*)filePriv->Get("T2");
   if (treePriv==0) {
     cout <<"Private tree not found!\n";
     return 0;
   }
   IOPrivateTree* privtree = new IOPrivateTree(treePriv);
   privtree->fill(ioc);

// eventwise tree
   printf("event tree ...\n");
   TFile* fileEvt = new TFile(sfileevt);
   IOEventTree* evttree = new IOEventTree(fileEvt);

// open ps output file
   gSystem->Exec("rm -f "+psname);
   TCanvas *can = new TCanvas("can", "Alignment");
   can->Print(psname+"[","ps"); 

   Analysis* ioc2 = new Analysis(ioc);

 //  --------------------------------------------------------------------------
 //  MY SETTINGS --------------------------------------------------------------
 //  --------------------------------------------------------------------------

  if (PxZoom) {
    ioc2->analysis(rootFile,can,psname,"PXB" ,floatVar,true,MinHit,ms,mr,itp,nbin,IterFile);
//    ioc2->analysis(rootFile,can,psname,"PXEC",floatVar,true,MinHit,ms,mr,itp,nbin,IterFile);
 //   ioc2->analysis(rootFile,can,psname,"PX"  ,floatVar,true,MinHit,ms,mr,itp,nbin,IterFile);
  }
  else {
    ioc2->analysis(rootFile,can,psname,"TIB_DS",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);    
    ioc2->analysis(rootFile,can,psname,"TIB_SS",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    ioc2->analysis(rootFile,can,psname,"TOB_DS",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    ioc2->analysis(rootFile,can,psname,"TOB_SS",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    // ioc2->analysis(rootFile,can,psname,"TIB",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    // ioc2->analysis(rootFile,can,psname,"TOB",floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    // ioc2->analysis(rootFile,can,psname,"TID"   ,floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    // ioc2->analysis(rootFile,can,psname,"TEC"   ,floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
    // ioc2->analysis(rootFile,can,psname,"STRIP" ,floatVar,false,MinHit,ms,mr,itp,nbin,IterFile);
  }

  evttree->analysis(rootFile,can,psname);


  can->Print(psname+"]","ps"); 
  rootFile->Write() ;
  rootFile->Close() ;
  printf("Finished.\n");
  return 0;
}



//-----------------------------------------------------------------------------

void PlotTLHisto(TH2F* histo, Option_t* opt)
{
  histo->SetStats(kFALSE);

  //   Int_t cullay[13]={0,9,15,18,22,24,27,30,32,36,39,45,54};
  Int_t cullay[13]={-27,-18,-12,-9,-5,-3,0,3,5,9,12,18,27};
     TString modname[12]={"TEC-","TOB-","TID-","TIB-","PBwd","PBar-",
                        "PBar+","PFwd","TIB+","TID+","TOB+","TEC+"};

   histo->Draw(opt);
   gPad->Update();
   Float_t ymax=gPad->GetUymax();
   Float_t ymin=gPad->GetUymin();
   //cout <<"ymin,max: " << ymin <<" "<<ymax << endl;
   for(Int_t k=0;k<12;k++){
     TLine *lin= new TLine(cullay[k],ymin,cullay[k],ymax);
     lin->SetLineWidth(1);
     lin->SetLineStyle(3);
     lin->SetLineColor(4);
     if (k>0) lin->Draw();
     //TLine lin(cullay[k],ymin,cullay[k],ymax);
     //lin.Draw();

     TLatex t;
     t.SetTextSize(0.04);
     t.SetTextAngle(90);
     //t.SetTextColor(10);
     Float_t tpos=cullay[k]+0.75*(cullay[k+1]-cullay[k]);
     t.DrawLatex(tpos,ymin+0.025*(ymax-ymin),modname[k]);
   }

}

//-----------------------------------------------------------------------------

void PlotTLHisto(TH1F* histo)
{
  histo->SetStats(kFALSE);

  //   Int_t cullay[13]={0,9,15,18,22,24,27,30,32,36,39,45,54};
  Int_t cullay[13]={-27,-18,-12,-9,-5,-3,0,3,5,9,12,18,27};
     TString modname[12]={"TEC-","TOB-","TID-","TIB-","PBwd","PBar-",
                        "PBar+","PFwd","TIB+","TID+","TOB+","TEC+"};

   histo->Draw("");
   gPad->Update();
   Float_t ymax=gPad->GetUymax();
   Float_t ymin=gPad->GetUymin();
   for(Int_t k=0;k<12;k++){
     TLine *lin= new TLine(cullay[k],ymin,cullay[k],ymax);
    lin->SetLineWidth(1);
     lin->SetLineStyle(3);
     lin->SetLineColor(4);
     if (k>0) lin->Draw();

     TLatex t;
     t.SetTextSize(0.04);
     t.SetTextAngle(90);
     Float_t tpos=cullay[k]+0.75*(cullay[k+1]-cullay[k]);
     t.DrawLatex(tpos,ymin+0.025*(ymax-ymin),modname[k]);
   }

}

