#ifndef IOAlignableDataTree_h
#define IOAlignableDataTree_h

#include <TROOT.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include "IOCombined.h"

class IOAlignableDataTree {

 public:

  // constructor
  IOAlignableDataTree(TTree* t);

  // destructor
  ~IOAlignableDataTree() {};

  // dump method
  void dump(void);
  // size
  Int_t size(void) { return (Int_t)fChain->GetEntries(); };
  // fill
  void fill(IOCombined* ioc, Int_t flag, Int_t iter);


 private:

  TVectorD eulerAngles(TMatrixD orig,Int_t flag);
  TMatrixD rot2mat(Double_t* rotmat);
  TMatrixD diffRot(TMatrixD rot1, TMatrixD rot2);
  TMatrixD toLocal(TMatrixD rot, TMatrixD absrot);
  TVectorD toLocal(TVectorD pos, TVectorD apos, TMatrixD arot);

  TTree* fChain;
  static Int_t Nent;

//Declaration of leaves types
   Int_t          Id;
   Double_t        Pos[3];
   Double_t        Rot[9];
   Int_t           ObjId;

//List of branches
   TBranch        *b_Id;   //!
   TBranch        *b_Pos;   //!
   TBranch        *b_Rot;   //!
   TBranch        *b_ObjId;   //!
};

// constructor ----------------------------------------------------------------

IOAlignableDataTree::IOAlignableDataTree(TTree* t)
{
  fChain = t;

  fChain->SetBranchAddress("Id",&Id);
  fChain->SetBranchAddress("ObjId",&ObjId);
  fChain->SetBranchAddress("Pos",Pos);
  fChain->SetBranchAddress("Rot",Rot);

  b_Id    = fChain->GetBranch("Id");
  b_ObjId = fChain->GetBranch("ObjId");
  b_Pos   = fChain->GetBranch("Pos");
  b_Rot   = fChain->GetBranch("Rot");
}

// dump -----------------------------------------------------------------------

void IOAlignableDataTree::dump(void)
{

  Int_t Nent=(Int_t)fChain->GetEntries();
  printf("\nIOAlignableDataTree has %d entries\n",Nent);
  printf("    Id   ObjId      pos[0]       pos[1]       pos[2] \n");

  for(Int_t Ient=0; Ient<Nent; Ient++) {
    fChain->GetEntry(Ient);
    printf("%8d %3d: %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f  \n",
	   Id,ObjId,Pos[0],Pos[1],Pos[2],Pos[3],Pos[4],Pos[5]);
  }
}

// fill -----------------------------------------------------------------------

void IOAlignableDataTree::fill(IOCombined* ioc, Int_t flag, Int_t iter)
{
  Int_t Nent=(Int_t)fChain->GetEntries();
  //printf("IOAlignableDataTree has %d entries\n",Nent);

//  cout << "[IOAlignableDataTree::fill()]\tNAli " <<ioc->NAli << endl ;
  for (Int_t i=0;i<ioc->NAli;i++) ioc->AliTest[i]=false;

  Int_t Iter=iter;

  for(Int_t Ient=0; Ient<Nent; Ient++) {
    fChain->GetEntry(Ient);

    int ind = ioc->findIndex(Id,ObjId);
    if (ind == -1 && flag==1) { ind=ioc->NAli; ioc->NAli++; }    

    // original position 
    if (flag ==1) { 
      ioc->AliId[ind]=Id;
      ioc->AliObjId[ind]=ObjId;

      // global
      ioc->AbsPosO[ind][0]=Pos[0];
      ioc->AbsPosO[ind][1]=Pos[1];
      ioc->AbsPosO[ind][2]=Pos[2];
      ioc->AbsPosO[ind][3]=0.; 
      ioc->AbsPosO[ind][4]=0.;
      ioc->AbsPosO[ind][5]=0.;
      // local
      ioc->AbsPosOl[ind][0]=0.;
      ioc->AbsPosOl[ind][1]=0.;
      ioc->AbsPosOl[ind][2]=0.;
      ioc->AbsPosOl[ind][3]=0.; 
      ioc->AbsPosOl[ind][4]=0.;
      ioc->AbsPosOl[ind][5]=0.;

      // store rotation matrix
      ioc->RotMatrixOArr()->AddAt(new TMatrixD(rot2mat(Rot)),ind);
      TVectorD p(3); p[0]=Pos[0]; p[1]=Pos[1]; p[2]=Pos[2];
      ioc->PosVectorOArr()->AddAt(new TVectorD(p),ind);
    }

    // misaligned positions
    else if (flag ==2 && ind>-1) { 
      TMatrixD* rotold2 = (TMatrixD*)(ioc->RotMatrixOArr()->At(ind));
      TMatrixD rotnew2 = rot2mat(Rot);
      TVectorD* posold2 = (TVectorD*)(ioc->PosVectorOArr()->At(ind));
      TVectorD posnew2(3); posnew2[0]=Pos[0]; posnew2[1]=Pos[1]; posnew2[2]=Pos[2];
      // global
      ioc->AbsPosM[ind][0]=Pos[0];
      ioc->AbsPosM[ind][1]=Pos[1];
      ioc->AbsPosM[ind][2]=Pos[2];
      TVectorD e2 = eulerAngles(diffRot(*rotold2,rotnew2),0);
      ioc->AbsPosM[ind][3]=e2[0];
      ioc->AbsPosM[ind][4]=e2[1];
      ioc->AbsPosM[ind][5]=e2[2];
      cout << posnew2[0] << " " << posnew2[1] << " " << posnew2[2] << endl;
      TVectorD myposold(*posold2);
      cout << myposold[0] << " " << myposold[1] << " " << myposold[2] << endl;
      // local
      TVectorD l2=toLocal(posnew2,*posold2,*rotold2)
	-toLocal(*posold2,*posold2,*rotold2);  // a che minchia serve???
      ioc->AbsPosMl[ind][0]=l2[0];
      ioc->AbsPosMl[ind][1]=l2[1];
      ioc->AbsPosMl[ind][2]=l2[2];
      TVectorD e2l=eulerAngles(toLocal(diffRot(*rotold2,rotnew2),*rotold2),0);
      ioc->AbsPosMl[ind][3]=e2l[0];
      ioc->AbsPosMl[ind][4]=e2l[1];
      ioc->AbsPosMl[ind][5]=e2l[2];
      // fill misaligned pos as iteration zero
      for (Int_t i=0; i<6;i++) { 
	ioc->AbsPosA[ind][0][i]  = ioc->AbsPosM[ind][i];
	ioc->AbsPosAl[ind][0][i] = ioc->AbsPosMl[ind][i];
      }
    }

    // aligned positions
    else if (flag ==3 && ind>-1 ) {
      TMatrixD* rotold = (TMatrixD*) (ioc->RotMatrixOArr()->At(ind));
      TMatrixD rotnew = rot2mat(Rot);
      TVectorD* posold = (TVectorD*)(ioc->PosVectorOArr()->At(ind));
      TVectorD posnew(3); posnew[0]=Pos[0]; posnew[1]=Pos[1]; posnew[2]=Pos[2];
     // global
      ioc->AbsPosA[ind][Iter][0]=Pos[0];
      ioc->AbsPosA[ind][Iter][1]=Pos[1];
      ioc->AbsPosA[ind][Iter][2]=Pos[2];
      TVectorD e = eulerAngles(diffRot(*rotold,rotnew),0);
      ioc->AbsPosA[ind][Iter][3]=e[0];
      ioc->AbsPosA[ind][Iter][4]=e[1];
      ioc->AbsPosA[ind][Iter][5]=e[2];
      //local
      TVectorD l=toLocal(posnew,*posold,*rotold)
	-toLocal(*posold,*posold,*rotold); // a che minchia serve???
      ioc->AbsPosAl[ind][Iter][0]=l[0];
      ioc->AbsPosAl[ind][Iter][1]=l[1];
      ioc->AbsPosAl[ind][Iter][2]=l[2];
      TVectorD el = eulerAngles(toLocal(diffRot(*rotold,rotnew),*rotold),0);
      ioc->AbsPosAl[ind][Iter][3]=el[0];
      ioc->AbsPosAl[ind][Iter][4]=el[1];
      ioc->AbsPosAl[ind][Iter][5]=el[2];

      ioc->AliUse[ind]=1;
      //cout << "Setting aliuse[" << ind <<"]" << endl ; 
      if (Iter>ioc->NIter) ioc->NIter=Iter;
      ioc->AliTest[ind]=true;
    }
  }


  // set those which habe not been aligned in this iteration
  // to values from previous iteration
  if (flag ==3) {
    for (Int_t i=0;i<ioc->NAli;i++) {
      if (ioc->AliUse[i]==true && ioc->AliTest[i]==false) {
	for (Int_t j=0;j<6;j++) 
          ioc->AbsPosA[i][Iter][j] = ioc->AbsPosA[i][Iter-1][j];
	printf("re-using values from prev. iter for %d \n",i);
      }
    }
  }

  //  printf("filled iteration %d \n",Iter);
  //printf("this iter, all Alignables, all Iterations: %d %d %d \n",
  //Iter,ioc->NAli,ioc->NIter);
  //printf("Iterations: %d \n",NIter);

}

//-----------------------------------------------------------------------------
// find R for which rot2 = rot1*R

TMatrixD IOAlignableDataTree::diffRot(TMatrixD rot1, TMatrixD rot2)
{
  Double_t det;
  TMatrixD result = rot1.Invert(&det) * rot2;
  return result;
}


//-----------------------------------------------------------------------------
// global -> local rotation matrix

TMatrixD IOAlignableDataTree::toLocal(TMatrixD rot, TMatrixD absrot)
{
  TMatrixD absrot2 = absrot;
  absrot2.Transpose(absrot);
  TMatrixD result = absrot * rot * absrot2;
  return result;
}

//-----------------------------------------------------------------------------
// global -> local position

TVectorD IOAlignableDataTree::toLocal(TVectorD pos, TVectorD apos, TMatrixD arot)
{
  TVectorD dpos(pos-apos);
  dpos *= arot;
  return dpos;
}

//-----------------------------------------------------------------------------

TMatrixD IOAlignableDataTree::rot2mat(Double_t* rotmat)
{
  TMatrixD  mat(3,3);

  mat[0][0]=rotmat[0];
  mat[0][1]=rotmat[1];
  mat[0][2]=rotmat[2];
  mat[1][0]=rotmat[3];
  mat[1][1]=rotmat[4];
  mat[1][2]=rotmat[5];
  mat[2][0]=rotmat[6];
  mat[2][1]=rotmat[7];
  mat[2][2]=rotmat[8];
  return mat;

}


//-----------------------------------------------------------------------------

TVectorD IOAlignableDataTree::eulerAngles(TMatrixD orig,Int_t flag)
{
  Double_t PI = 3.1415927;

  //  AlgebraicMatrix orig =  algebraicMatrix(rot);
  TMatrixD testangle(3,2);
  TVectorD returnangle(3);

  returnangle[0]=0.;
  returnangle[1]=0.;
  returnangle[2]=0.;

  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<3; j++) {
      if (fabs(orig[i][j])>1.001)  { 
        cout <<"Euler: Error in RotMat! i,j,val: " << i <<","<<j<<","<<orig[i][j]<<endl;
        return(returnangle);
      }
    }
  }
 
  if(orig[2][0]!=1.0) // If angle1 is not +-PI/2
    {

      if(flag==0) // assuming -PI/2 < angle1 < PI/2 
	{
//	  testangle[1][flag]=asin(-orig[2][0]);
	  testangle[1][flag]=asin(orig[2][0]);
	}

      if(flag==1) // assuming angle1 < -PI/2 or angle1 >PI/2
	{
//	  testangle[1][flag]=PI-asin(-orig[2][0]);
	  testangle[1][flag]=PI-asin(orig[2][0]);
	}


      if(cos(testangle[1][flag])*orig[2][2]>0)
	{
	  testangle[0][flag]=atan(-orig[2][1]/orig[2][2]);
	}
      else
	{
	  testangle[0][flag]=atan(-orig[2][1]/orig[2][2])+PI;
	}

      if(cos(testangle[1][flag])*orig[0][0]>0)
	{
	  testangle[2][flag]=atan(-orig[1][0]/orig[0][0]);
	}
      else
	{
	  testangle[2][flag]=atan(-orig[1][0]/orig[0][0])+PI;
	}
    }

  else // if angle1 == +-PI/2
    {
      testangle[1][flag]=PI/2; // chose positve Solution 
      if(orig[2][2]>0)
	{
        testangle[2][flag]=atan(orig[1][2]/orig[1][1]);
        testangle[0][flag]=0;
	  }
    }
 
  for(int i=0;i<3;i++)
    {
      returnangle[i]=testangle[i][flag];
    }
  //  cout << returnangle[0]<<","<< returnangle[1]<<","<<returnangle[2]<<","<<endl;
  return(returnangle);
}

//-----------------------------------------------------------------------------



#endif
