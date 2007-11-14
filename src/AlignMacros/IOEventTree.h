#ifndef IOEventTree_h
#define IOEventTree_h

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>

//-----------------------------------------------------------------------------

class IOEventTree {

 public:

   IOEventTree(TFile* file); // constructor
   ~IOEventTree() {}; // destructor
   Int_t initTree(int iter);
   void analysis(TFile * rootFile, TCanvas* can, TString psname); // fill

 private:

   TTree* fChain;
   TFile* theFile;

   static const Int_t TMAX=99;

// Declaration of leaves types
   Int_t           Run;
   Int_t           Event;
   Int_t           Ntracks;
   Int_t           Nhits[TMAX];   //[Ntracks]
   Float_t         Pt[TMAX];   //[Ntracks]
   Float_t         PtTPE[TMAX];   //[Ntracks]
   Float_t         Eta[TMAX];   //[Ntracks]
   Float_t         Phi[TMAX];   //[Ntracks]
   Float_t         Theta[TMAX];   //[Ntracks]
   Float_t         Cov[TMAX][6];   //[Ntracks]
   Float_t         M12;
//   Float_t         M12_init;
   Double_t        M12_estim;
   Float_t         Chi2[TMAX];
   Float_t         Chi2n[TMAX];
   Int_t           Dof[TMAX];

// List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_Nhits;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_PtTPE;   //!
   TBranch        *b_Eta;   //!
   TBranch        *b_Phi;   //!
   TBranch        *b_Theta;   //!
   TBranch        *b_Cov;   //!
   TBranch        *b_M12;   //!
//   TBranch        *b_M12_init;   //!
   TBranch        *b_M12_estim;   //!
   TBranch        *b_Chi2;   //!
   TBranch        *b_Chi2n;   //!
   TBranch        *b_Dof;   //!
  
};

// constructor ----------------------------------------------------------------

IOEventTree::IOEventTree(TFile* file)
{
   theFile = file;
   //Int_t iret = initTree(1);
   //if (iret==-1) cout <<"ERROR: tree not found for iter \n";
}

Int_t IOEventTree::initTree(int iter)
{

   TString tname="T1";
   char iterString[5];
   sprintf(iterString, "%i",iter);
   tname.Append("_");
   tname.Append(iterString);

   TTree* tree = (TTree*)theFile->Get(tname);

   if (tree==0) return -1;

   fChain = tree;

//   fChain->SetBranchAddress("Run",&Run);
//   fChain->SetBranchAddress("Event",&Event);
   fChain->SetBranchAddress("Ntracks",&Ntracks);
   fChain->SetBranchAddress("Nhits",Nhits);
   fChain->SetBranchAddress("Pt",Pt);
//   fChain->SetBranchAddress("PtTPE",PtTPE);
   fChain->SetBranchAddress("Eta",Eta);
   fChain->SetBranchAddress("Phi",Phi);
//   fChain->SetBranchAddress("Theta",Theta);
//   fChain->SetBranchAddress("Cov",Cov);
//   fChain->SetBranchAddress("M12",&M12);
//   fChain->SetBranchAddress("M12_init",&M12_init);
//   fChain->SetBranchAddress("M12_estim",&M12_estim);
//   fChain->SetBranchAddress("Chi2",Chi2);
   fChain->SetBranchAddress("Chi2n",Chi2n);
//   fChain->SetBranchAddress("Dof",Dof);

//   b_Run = fChain->GetBranch("Run");
//   b_Event = fChain->GetBranch("Event");
   b_Ntracks = fChain->GetBranch("Ntracks");
   b_Nhits = fChain->GetBranch("Nhits");
   b_Pt = fChain->GetBranch("Pt");
//   b_PtTPE = fChain->GetBranch("PtTPE");
   b_Eta = fChain->GetBranch("Eta");
   b_Phi = fChain->GetBranch("Phi");
//   b_Theta = fChain->GetBranch("Theta");
//   b_Cov = fChain->GetBranch("Cov");
//   b_M12 = fChain->GetBranch("M12");
//   b_M12_init = fChain->GetBranch("M12_init");
//   b_M12_estim = fChain->GetBranch("M12_estim");
//   b_Chi2 = fChain->GetBranch("Chi2");
   b_Chi2n = fChain->GetBranch("Chi2n");
//   b_Dof = fChain->GetBranch("Dof");

   return 0;
}

// fill -----------------------------------------------------------------------

void IOEventTree::analysis(TFile * rootFile, TCanvas* can, TString psname)
{

   //  static const Int_t MAXITER=41;
   //  static const Int_t MAXITER=51;
   static const Int_t MAXITER=98;

  // book histos --------------------------------------------------------------

  TH1F* hntracks = new TH1F("hntracks"," N(Tracks) ; N(Tracks)",10,1,10);
  TH1F* hnh      = new TH1F("hnh"," Hits(Track) ; Hits(Track)",60,0,60);
  TH1F* hpt      = new TH1F("hpt"," P_{T}(Track) ; P_{T} / GeV",50,0,200);
//  TH1F* hpt      = new TH1F("hpt"," P_{T}(Track) ; P_{T} / GeV",50,0,10);
  TH1F* hptTPE   = new TH1F("hptTPE"," P_{T}TPE(Track) ; P_{T}TPE / GeV",50,0,200);
  TH1F* heta     = new TH1F("heta"," #eta(Track) ; #eta",50,-3.0,3.0);
//  TH1F* hM12     = new TH1F("hM12"," M12 ",80,70.,110.);
//  TH1F* hM12_init = new TH1F("hM12_init"," M12_init ",80,70.,110.);
//  TH1F* hM12_estim = new TH1F("hM12_estim"," M12_estim ",80,70.,110.);

  TH1F hchi2[MAXITER];
  TH1F hM12[MAXITER];
  TH1F hM12estim[MAXITER];
  

  char iterString[2];
  for(Int_t i=1; i<MAXITER; i++) {

    sprintf(iterString, "%d",i);
    
    TString tname="Chi2 for Iteration ";
    tname.Append(iterString);
    tname.Append(" ; Chi2(Track)");
    TString histoName = "hchi2_" ;
    histoName.Append(iterString) ;
    
    /* TString tname2="m12 for Iteration ";
    tname2.Append(iterString);
    tname2.Append(" ; InvMass");
    TString histoName2 = "hM12_" ;
    histoName2.Append(iterString) ;
   
    TString tname3="m12 refit for Iteration ";
    tname3.Append(iterString);
    tname3.Append(" ; InvMass");
    TString histoName3 = "hM12estim_" ;
    histoName3.Append(iterString) ;  */

    hchi2[i]=TH1F(histoName,tname,50,0.,5.);
    // hM12[i] =TH1F(histoName2, tname2, 80, 70.,110.) ;
    // hM12estim[i] =TH1F(histoName3, tname3, 80, 70.,110.) ;
    
  }

  // fill histos --------------------------------------------------------------

  Int_t iret,Nent,itermin;

  itermin=1;
  iret = initTree(1);
  if (iret==-1) { 
    cout <<"ERROR: event tree not found for iter1! try iter2 ...\n";
    iret = initTree(2);
    if (iret==-1) { 
      cout <<"ERROR: event tree not found for iter2! EXIT! \n";
      return;
    }
    itermin=2;
  }

  Nent=(Int_t)fChain->GetEntries();
  printf("IOEventTree has %d entries\n",Nent);

  int ichecktpe = 0;
  int ichecktot = 0;
  for(Int_t Ient=0; Ient<Nent; Ient++) {
    fChain->GetEntry(Ient);
    hntracks->Fill(Ntracks);
    for(Int_t Itr=0; Itr<Ntracks; Itr++) {
      hnh->Fill(Nhits[Itr]);
      hpt->Fill(Pt[Itr]);
      ichecktot++;
      heta->Fill(Eta[Itr]);
      if (PtTPE[Itr] != 0) {
	hptTPE->Fill(PtTPE[Itr]);
	ichecktpe++;
      }
    }
  }
  cout << " tot " << ichecktot << endl;
  // cout << " tpe " << ichecktpe << endl;

  
  Int_t maxiter=0;
  for(Int_t i=itermin; i<MAXITER; i++) {
    iret = initTree(i);
    if (iret==-1) { maxiter=i-1; break; }
    for(Int_t Ient=0; Ient<(fChain->GetEntries()); Ient++) {
      fChain->GetEntry(Ient);
      for(Int_t Itr=0; Itr<Ntracks; Itr++) hchi2[i].Fill(Chi2n[Itr]);
      // hM12[i].Fill(M12);
      // if (M12_estim != 0) hM12estim[i].Fill(M12_estim);
    }
    if(i == MAXITER-1 ) {maxiter = i ; break ;}
    cout << "done" << endl ; 
  }

  cout <<"Iterations filled: " << maxiter << endl;

  // make plots ---------------------------------------------------------------


   can->Clear();
   can->Divide(2,2);

   can->cd(1); hntracks->Draw();
   can->cd(2); hnh->Draw();
   can->cd(3); hpt->Draw();
   can->cd(4); hptTPE->Draw();
//   can->cd(4); heta->Draw();
//   can->cd(4); hM12->Draw();

   can->Update();
   can->Print(psname,"ps"); 

   // - - - - - - - - - - - - - - - - -

//   can->Clear();
//   can->Divide(2,2);
//   can->cd(1); hM12->Draw();
//   can->cd(2); hM12_init->Draw();
//   can->cd(2); hM12_estim->Draw();

//   can->Update();
//   can->Print(psname,"ps"); 

   // - - - - - - - - - - - - - - - - -
  //chi2
   
   can->Clear();
   can->Divide(3,3);
   for(Int_t i=1; i<=9 && i<=maxiter; i++) { 
     can->cd(i);  hchi2[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 

   if (maxiter>9) {
   can->Clear();
   can->Divide(3,3);
   for(Int_t i=10; i<=18 && i<=maxiter; i++) { 
     can->cd(i-9);  hchi2[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   /* if (maxiter>18) {
   can->Clear();
   can->Divide(3,3);
   for(Int_t i=19; i<=27 && i<=maxiter; i++) { 
     can->cd(i-18);  hchi2[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>27) {
   can->Clear();
   can->Divide(3,3);
   for(Int_t i=28; i<=36 && i<=maxiter; i++) { 
     can->cd(i-27);  hchi2[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   } */

   // - - - - - - - - - - - - - - - - -
  //m12
   
   /* can->Clear();
   can->Divide(3,2);
   for(Int_t i=1; i<=3 && i<=maxiter; i++) { 
     can->cd(i);  hM12[i].Draw();
     can->cd(i+3);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 

   if (maxiter>3) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=4; i<=6 && i<=maxiter; i++) { 
     can->cd(i-3);  hM12[i].Draw();
     can->cd(i);    hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>6) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=7; i<=9 && i<=maxiter; i++) { 
     can->cd(i-6);  hM12[i].Draw();
     can->cd(i-3);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>9) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=10; i<=12 && i<=maxiter; i++) { 
     can->cd(i-9);  hM12[i].Draw();
     can->cd(i-6);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>12) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=13; i<=15 && i<=maxiter; i++) { 
     can->cd(i-12);  hM12[i].Draw();
     can->cd(i-9);   hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>15) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=16; i<=18 && i<=maxiter; i++) { 
     can->cd(i-15);  hM12[i].Draw();
     can->cd(i-12);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>18) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=19; i<=21 && i<=maxiter; i++) { 
     can->cd(i-18);  hM12[i].Draw();
     can->cd(i-15);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }

   if (maxiter>21) {
   can->Clear();
   can->Divide(3,2);
   for(Int_t i=22; i<=24 && i<=maxiter; i++) { 
     can->cd(i-21);  hM12[i].Draw();
     can->cd(i-18);  hM12estim[i].Draw();
   }
   can->Update();
   can->Print(psname,"ps"); 
   }*/

}

//-----------------------------------------------------------------------------





#endif
