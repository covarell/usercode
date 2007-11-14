#ifndef IOCombined_h
#define IOCombined_h

class IOCombined {

 public:

  // constructor
  IOCombined() : NAli(0), NIter(0) 
  {
    for (Int_t i=0;i<ALIMAX;i++) AliUse[i]=0;
    cout <<"creating pos/rot array of size " << ALIMAX << endl;
    theRotMatrixOArr = new TObjArray(ALIMAX,0);
    thePosVectorOArr = new TObjArray(ALIMAX,0);

  };

  // destructor
  ~IOCombined() {delete theRotMatrixOArr; };

  TObjArray* RotMatrixOArr(void) {return theRotMatrixOArr; };
  TObjArray* PosVectorOArr(void) {return thePosVectorOArr; };

  // helper methods
  Int_t findIndex(Int_t id, Int_t objid);
  Bool_t isDoubleSided(Int_t iali);
  Int_t getlayer(Int_t iali);

  // data members

  static const Int_t ALIMAX=14000;
//  static const Int_t ITERMAX=41;
//  static const Int_t ITERMAX=51;
  static const Int_t ITERMAX=98;

  Int_t   AliUse[ALIMAX];
  Int_t   AliId[ALIMAX];
  Int_t   AliObjId[ALIMAX];
  Double_t AbsPosO[ALIMAX][6];
  Double_t AbsPosOl[ALIMAX][6];
  Double_t AbsPosM[ALIMAX][6];
  Double_t AbsPosMl[ALIMAX][6];
  Double_t AbsPosA[ALIMAX][ITERMAX][6];
  Double_t AbsPosAl[ALIMAX][ITERMAX][6];
  Int_t    AliNhit[ALIMAX];
  Double_t AliPosR[ALIMAX],AliPosZ[ALIMAX],AliPosPhi[ALIMAX];
  Int_t    AliType[ALIMAX],AliLayer[ALIMAX];
  Bool_t   AliTest[ALIMAX];
  Double_t AliPar[ALIMAX][ITERMAX][6];
  Double_t AliParErr[ALIMAX][ITERMAX][6];
  TObjArray* theRotMatrixOArr;
  TObjArray* thePosVectorOArr;

  Int_t NAli;
  Int_t NIter;

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

Int_t IOCombined::findIndex(Int_t id, Int_t objid)
{
  for (Int_t i=0; i<NAli; i++) 
    if (AliId[i] == id && AliObjId[i] == objid) return i;
  return -1;
}

//-----------------------------------------------------------------------------

Bool_t IOCombined::isDoubleSided(Int_t iali)
{
  Int_t imod=AliType[iali];
  Int_t ilay=AliLayer[iali];
  // Barrel: TIB+TOB
  if (TMath::Abs(imod)==3 || TMath::Abs(imod)==5) if (ilay<=2) return true;
  // pixel barrel or endcap
  if (TMath::Abs(imod)==1) return true;
  if (TMath::Abs(imod)==2) return true;

  return false;
}

//-----------------------------------------------------------------------------

Int_t IOCombined::getlayer(Int_t iali)
{

  Int_t imod=AliType[iali];
  Int_t ilay=AliLayer[iali];
  Float_t theZ = AbsPosO[iali][2];
 
  Int_t cullay[6]={0,3,5,9,12,18};

  if (imod<0 || imod>6 || ilay>99 || ilay<0 ) return(-1);
  
  Int_t ilayer;
  if (theZ > 0) {
    ilayer = cullay[imod-1]+ilay-1;
  } else {
    ilayer = -cullay[imod-1]-ilay;
  }

  return(ilayer);

}




#endif



