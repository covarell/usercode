#ifndef IOPrivateTree_h
#define IOPrivateTree_h

#include <TROOT.h>
#include <TTree.h>

//#include "IOAlignmentTree.h"
#include "IOCombined.h"

class IOPrivateTree {

 public:

  // constructor
  IOPrivateTree(TTree* t);

  // destructor
  ~IOPrivateTree();

  // size
  Int_t size(void);
  // fill
  void fill(IOCombined* ioc);

 private:

  TTree* fChain;
  static Int_t Nent;

//Declaration of leaves types
   Int_t           Nhit;
   Int_t           Niter;
   Int_t           Type;
   Int_t           Layer;
   Float_t         Xpos;
   Float_t         Ypos;
   Float_t         Zpos;
   Float_t         Eta;
   Float_t         Phi;
   Int_t           Id;
   Int_t           ObjId;
   Double_t        AbsPosO[3];
   Double_t        AbsRotO[9];
   Double_t        AbsPosM[3];
   Double_t        AbsRotM[9];
   Double_t        RelPosM[3];
   Double_t        RelRotM[9];
   Double_t        AbsPosA[3];
   Double_t        AbsRotA[9];
   Double_t        RelPosA[3];
   Double_t        RelRotA[9];
   Double_t        AbsPosI[99][3];   //[Niter]
   Double_t        AbsRotI[99][9];   //[Niter]
   Double_t        LocParI[99][6];   //[Niter]


//List of branches
   TBranch        *b_Nhit;   //!
   TBranch        *b_Niter;   //!
   TBranch        *b_Type;   //!
   TBranch        *b_Layer;   //!
   TBranch        *b_Xpos;   //!
   TBranch        *b_Ypos;   //!
   TBranch        *b_Zpos;   //!
   TBranch        *b_Eta;   //!
   TBranch        *b_Phi;   //!
   TBranch        *b_Id;   //!
   TBranch        *b_ObjId;   //!
   TBranch        *b_AbsPosO;   //!
   TBranch        *b_AbsRotO;   //!
   TBranch        *b_AbsPosM;   //!
   TBranch        *b_AbsRotM;   //!
   TBranch        *b_RelPosM;   //!
   TBranch        *b_RelRotM;   //!
   TBranch        *b_AbsPosA;   //!
   TBranch        *b_AbsRotA;   //!
   TBranch        *b_RelPosA;   //!
   TBranch        *b_RelRotA;   //!
   TBranch        *b_AbsPosI;   //!
   TBranch        *b_AbsRotI;   //!
   TBranch        *b_LocParI;   //!

};

// constructor ----------------------------------------------------------------

IOPrivateTree::IOPrivateTree(TTree* t)
{
  fChain = t;

   fChain->SetBranchAddress("Nhit",&Nhit);
//   fChain->SetBranchAddress("Niter",&Niter);
   fChain->SetBranchAddress("Type",&Type);
   fChain->SetBranchAddress("Layer",&Layer);
   fChain->SetBranchAddress("Xpos",&Xpos);
   fChain->SetBranchAddress("Ypos",&Ypos);
   fChain->SetBranchAddress("Zpos",&Zpos);
   fChain->SetBranchAddress("Eta",&Eta);
   fChain->SetBranchAddress("Phi",&Phi);
   fChain->SetBranchAddress("Id",&Id);
   fChain->SetBranchAddress("ObjId",&ObjId);
//   fChain->SetBranchAddress("AbsPosO",AbsPosO);
//   fChain->SetBranchAddress("AbsRotO",AbsRotO);
//   fChain->SetBranchAddress("AbsPosM",AbsPosM);
//   fChain->SetBranchAddress("AbsRotM",AbsRotM);
//   fChain->SetBranchAddress("RelPosM",RelPosM);
//   fChain->SetBranchAddress("RelRotM",RelRotM);
//   fChain->SetBranchAddress("AbsPosA",AbsPosA);
//   fChain->SetBranchAddress("AbsRotA",AbsRotA);
//   fChain->SetBranchAddress("RelPosA",RelPosA);
//   fChain->SetBranchAddress("RelRotA",RelRotA);
//   fChain->SetBranchAddress("AbsPosI",&AbsPosI);
//   fChain->SetBranchAddress("AbsRotI",&AbsRotI);
//   fChain->SetBranchAddress("LocParI",&LocParI);

   b_Nhit  = fChain->GetBranch("Nhit");
//   b_Niter = fChain->GetBranch("Niter");
   b_Type  = fChain->GetBranch("Type");
   b_Layer = fChain->GetBranch("Layer");
   b_Xpos  = fChain->GetBranch("Xpos");
   b_Ypos  = fChain->GetBranch("Ypos");
   b_Zpos  = fChain->GetBranch("Zpos");
   b_Eta   = fChain->GetBranch("Eta");
   b_Phi   = fChain->GetBranch("Phi");
   b_Id    = fChain->GetBranch("Id");
   b_ObjId = fChain->GetBranch("ObjId");
//   b_AbsPosO = fChain->GetBranch("AbsPosO");
//   b_AbsRotO = fChain->GetBranch("AbsRotO");
//   b_AbsPosM = fChain->GetBranch("AbsPosM");
//   b_AbsRotM = fChain->GetBranch("AbsRotM");
//   b_RelPosM = fChain->GetBranch("RelPosM");
//   b_RelRotM = fChain->GetBranch("RelRotM");
//   b_AbsPosA = fChain->GetBranch("AbsPosA");
//   b_AbsRotA = fChain->GetBranch("AbsRotA");
//   b_RelPosA = fChain->GetBranch("RelPosA");
//   b_RelRotA = fChain->GetBranch("RelRotA");
//   b_AbsPosI = fChain->GetBranch("AbsPosI");
//   b_AbsRotI = fChain->GetBranch("AbsRotI");
//   b_LocParI = fChain->GetBranch("LocParI");


   //  Int_t Nent=(Int_t)fChain->GetEntries();
  //  printf("IOPrivateTree created with %d entries.\n",Nent);

}

// destructor -----------------------------------------------------------------

IOPrivateTree::~IOPrivateTree()
{}

// size -----------------------------------------------------------------------

Int_t IOPrivateTree::size(void)
{
  return (Int_t)fChain->GetEntries();
}

// fill -----------------------------------------------------------------------

void IOPrivateTree::fill(IOCombined* ioc)
{
  Int_t Nent=(Int_t)fChain->GetEntries();
  printf("IOPrivateTree has %d entries\n",Nent);

  for(Int_t Ient=0; Ient<Nent; Ient++) {
    fChain->GetEntry(Ient);

    int ind = ioc->findIndex(Id,ObjId);
    if (ind == -1) {
       ind=ioc->NAli; ioc->NAli++; 
       printf("Warning: new alignable created ...\n");
    }    

    ioc->AliNhit[ind]=Nhit;

    ioc->AliPosR[ind]  =sqrt(Xpos*Xpos+Ypos*Ypos);
    ioc->AliPosZ[ind]  =Zpos;
    ioc->AliPosPhi[ind]=Phi;
    ioc->AliType[ind]=Type;
    ioc->AliLayer[ind]=Layer;

    // printf("filling hits: %d\n",AliNhit[ind]);
    //    printf("r,z,phi: %f %f %f\n",AliPosR[ind],AliPosZ[ind],AliPosPhi[ind]);
//    printf("type,layer: %d %d\n",Type,Layer);

  }
}

//-----------------------------------------------------------------------------





#endif
