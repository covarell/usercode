#ifndef IOAlignmentParameterTree_h
#define IOAlignmentParameterTree_h

#include <TROOT.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include "IOCombined.h"

class IOAlignmentParameterTree {

 public:

  // constructor
  IOAlignmentParameterTree(TTree* t);

  // destructor
  ~IOAlignmentParameterTree() {};

  // dump method
  void dump(void);
  // size
  Int_t size(void) { return (Int_t)fChain->GetEntries(); };
  // fill
  void fill(IOCombined* ioc, Int_t iter);


 private:

  TTree* fChain;
  static Int_t Nent;

// Declaration of leaves types
   Int_t           parSize;
   Int_t           Id;
   Double_t        Par[6];   //[CovRang]
   Int_t           covarSize;
   Double_t        Cov[21];   //[covarRang]
   Int_t           ObjId;

// List of branches
   TBranch        *b_CovRang;   //!
   TBranch        *b_Id;   //!
   TBranch        *b_Par;   //!
   TBranch        *b_covarRang;   //!
   TBranch        *b_Cov;   //!
   TBranch        *b_ObjId;   //!

};

// constructor ----------------------------------------------------------------

IOAlignmentParameterTree::IOAlignmentParameterTree(TTree* t)
{
   fChain = t;

   fChain->SetBranchAddress("parSize",&parSize);
   fChain->SetBranchAddress("Id",&Id);
   fChain->SetBranchAddress("Par",Par);
   fChain->SetBranchAddress("covarSize",&covarSize);
   fChain->SetBranchAddress("Cov",Cov);
   fChain->SetBranchAddress("ObjId",&ObjId);

   b_CovRang   = fChain->GetBranch("parSize");
   b_Id        = fChain->GetBranch("Id");
   b_Par       = fChain->GetBranch("Par");
   b_covarRang = fChain->GetBranch("covarSize");
   b_Cov       = fChain->GetBranch("Cov");
   b_ObjId     = fChain->GetBranch("ObjId");
}

// fill -----------------------------------------------------------------------

void IOAlignmentParameterTree::fill(IOCombined* ioc, Int_t iter)
{
  Int_t Nent=(Int_t)fChain->GetEntries();
  if (iter==1) printf("IOAlignmentParameterTree has %d entries\n",Nent);

  for(Int_t Ient=0; Ient<Nent; Ient++) {
    fChain->GetEntry(Ient);
    Int_t ind = ioc->findIndex(Id,ObjId);
    if (ind>-1) {
      for (Int_t i=0;i<6;i++) ioc->AliPar[ind][iter][i]=Par[i];
      ioc->AliParErr[ind][iter][0]=sqrt(Cov[ 0]);
      ioc->AliParErr[ind][iter][1]=sqrt(Cov[ 6]);
      ioc->AliParErr[ind][iter][2]=sqrt(Cov[11]);
      ioc->AliParErr[ind][iter][3]=sqrt(Cov[15]);
      ioc->AliParErr[ind][iter][4]=sqrt(Cov[18]);
      ioc->AliParErr[ind][iter][5]=sqrt(Cov[20]);
    }
  }
}

#endif

