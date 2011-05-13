#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "TH2.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"

#include "RooFit.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"

#include <TCanvas.h>
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TVector3.h>
#include <math.h>

using namespace RooFit;
using namespace std;

bool Cowboy(int mu1_charge, TLorentzVector* mu1, TLorentzVector *mu2){
  return (mu1_charge * mu1->DeltaPhi(*mu2) >0.);
} 

bool isAccept(const TLorentzVector* aMuon) {
   // use *OLD* muon kinematical cuts (eta dependent momentum / pT cuts )
   // return (fabs(aMuon->Eta()) < 2.4 &&
   //        ((fabs(aMuon->Eta()) < 1.3 && aMuon->Pt() > 3.3) ||
   //        (fabs(aMuon->Eta()) > 1.3 && fabs(aMuon->Eta()) < 2.2 && aMuon->P() > 2.9) ||
   //        (fabs(aMuon->Eta()) > 2.2 && aMuon->Pt() > 0.8)));

   // use *NEW* muon kinematical cuts (eta dependent momentum / pT cuts )
   return (fabs(aMuon->Eta()) < 2.4 &&
           ((fabs(aMuon->Eta()) < 1.2 && aMuon->Pt() > 4.0) ||
	    (fabs(aMuon->Eta()) > 1.2 && aMuon->Pt() > 3.3)));

   // *REMOVE* muon kinematical cuts (eta dependent momentum / pT cuts )
   // by just returning TRUE
   //  return true;
}

double CorrectMass(const TLorentzVector* mu1,const TLorentzVector* mu2, int mode){  
  double CMass=0;
  const double mumass=0.105658;
  double k1,k2;
  double pt1=mu1->Pt();
  double pt2=mu2->Pt();
  double eta1=mu1->Eta();
  double eta2=mu2->Eta();
  if (mode==1){
    k1=1.0009;//constant scale correction
    k2=1.0009;
  }
  if (mode==2){
    k1=1.0019-0.0004*pt1;
    k2=1.0019-0.0004*pt2; // pt dependent correction
  }
  if (mode==3){
      double a0=0.00038; //3.8 * pow(10,-4);
      double a1=0.0;
      double a2=0.0003; //3.0 * pow(10,-4);
      double a3=0.0;

      k1=1+a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=1+a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  if (mode == 4){
      double a0=1.002;
      double a1=-0.002;
      double a2=0.001;
      double a3=-0.0001;

      k1=a0+a1*fabs(eta1)+a2*eta1*eta1+a3*pt1;
      k2=a0+a1*fabs(eta2)+a2*eta2*eta2+a3*pt2;// pt and eta dependent
  }

  TVector3 mom1=mu1->Vect();
  TVector3 mom2=mu2->Vect();
  mom1=k1*mom1; 
  mom2=k2*mom2;
  double E1=sqrt(mom1.Mag2()+(mumass*mumass));
  double E2=sqrt(mom2.Mag2()+(mumass*mumass));
  TVector3 momtot=mom1+mom2;
  CMass=sqrt((E1+E2)*(E1+E2)-momtot.Mag2());
  return CMass;
}



int main(int argc, char* argv[]) {
  
  double JpsiMassMin = 2.5;
  double JpsiMassMax = 4.2;
  const double JpsiPtMin = 0;
  const double JpsiPtMax = 100;
  const double JpsiYMin =0;
  const double JpsiYMax = 2.4;
  const double JpsiCtMin = -2.0;
  const double JpsiCtMax = 3.5;
  
  char fileName[100];
  int doMerge = 0;
  int doWideRange = 0;
  
  if ( argc < 4 ){
    cout << "missing arguments: insert inputFile and mergePsiPrime" << endl; 
    return 1;
  }
  strcpy(fileName,argv[1]);
  doMerge = atoi(argv[2]);
  doWideRange = atoi(argv[3]);

  TFile *file= TFile::Open(fileName);
  TTree * Tree=(TTree*)file->Get("data");

  TLorentzVector* JP= new TLorentzVector;
  TLorentzVector* m1P= new TLorentzVector;
  TLorentzVector* m2P= new TLorentzVector;
  
  double vprob, theCt, theCtErr;
  int trig,trig2,theCat,Jq;

  static const unsigned int rapRegions = 3;
  float rapLimits[rapRegions+1] = {JpsiYMin,1.2,1.6,JpsiYMax};
  RooDataSet* dataJpsi[rapRegions];
  RooDataSet* dataPsip[rapRegions];
  RooRealVar* Jpsi_Mass;
  RooRealVar* Psip_Mass;      
  RooRealVar* Jpsi_Pt;
  RooRealVar* Jpsi_Ct;
  RooRealVar* Jpsi_CtErr;
  RooRealVar* Jpsi_Y;
  RooCategory* Jpsi_Type;
  RooCategory* Jpsi_PsiP;

  if (doWideRange) {
    JpsiMassMax += 0.5;
    Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  } else {
    Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",JpsiMassMin,3.5,"GeV/c^{2}");
  }
  if (doMerge) {
    Psip_Mass = new RooRealVar("PsiP_Mass","psi' mass",3.3,JpsiMassMax,"GeV/c^{2}");
  } else {
    Psip_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",3.3,JpsiMassMax,"GeV/c^{2}");
  }
  Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
  Jpsi_Y = new RooRealVar("Jpsi_Y","J/psi y",-JpsiYMax,JpsiYMax);
  Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi");
  Jpsi_PsiP = new RooCategory("Jpsi_PsiP","Jpsi or psiPrime");
  Jpsi_Ct = new RooRealVar("Jpsi_Ct","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
  Jpsi_CtErr = new RooRealVar("Jpsi_CtErr","J/psi ctau error",0.,4.,"mm");

  Jpsi_Type->defineType("GG",0);
  Jpsi_Type->defineType("GT",1);
  Jpsi_Type->defineType("TT",2);

  Jpsi_PsiP->defineType("J",0);
  Jpsi_PsiP->defineType("P",1);

  Tree->SetBranchAddress("JpsiP", &JP);
  Tree->SetBranchAddress("muPosP", &m1P);
  Tree->SetBranchAddress("muNegP", &m2P);
  Tree->SetBranchAddress("JpsiVprob", &vprob);
  Tree->SetBranchAddress("HLT_DoubleMu0", &trig);
  Tree->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &trig2);
  Tree->SetBranchAddress("JpsiType", &theCat);
  Tree->SetBranchAddress("JpsiCharge", &Jq);
  Tree->SetBranchAddress("Jpsict", &theCt);
  Tree->SetBranchAddress("JpsictErr", &theCtErr);

  RooArgList varlist(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr);
  RooArgList varlist2(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr);
  RooArgList varlist3(*Jpsi_Mass,*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_PsiP);
    
  for (unsigned int i = 0; i < rapRegions; i++) {
    if (doMerge) {
      dataJpsi[i] = new RooDataSet("dataAll","A sample",varlist3);
    } else {
      dataJpsi[i] = new RooDataSet("dataJpsi","A sample",varlist);
      dataPsip[i] = new RooDataSet("dataPsip","A sample",varlist2);
    }
  }

  int n = Tree->GetEntries();
  int nev=n;

  for (int ev=0; ev<nev; ++ev) {
    if (ev%5000==0) cout << ">>>>> EVENT " << ev << " / " << n<<  endl;
    
    Tree->GetEntry(ev);

    double theMass =JP->M();
    double CMass=CorrectMass(m1P,m2P,4);
    //cout << " Mass " << theMass << "   Corrected " << CMass << endl; 
    if (CMass!=0) theMass=CMass;

    double theRapidity=JP->Rapidity();
    double thePt=JP->Pt();
    /*
    cout << JP->Pt()<< endl;
    cout << JP->M()<< endl;
    cout << JP->Rapidity() << endl;
    // cout << "Triggered " << trig << endl;
    */
    bool ok1=isAccept(m1P);
    bool ok2=isAccept(m2P);
    
    // COWBOY CUTS!
    if (Cowboy(1,m1P,m2P) && fabs(JP->Rapidity()) > 1.6) continue;
    // if (Cowboy(1,m1P,m2P)) continue;
    // if (!Cowboy(1,m1P,m2P)) continue;

    //cout << "Accepted1 " << ok1 << "  Accepted2 " << ok2 << endl;
    if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
	theCt > JpsiCtMin && theCt < JpsiCtMax && 
	thePt > JpsiPtMin && thePt < JpsiPtMax && 
	fabs(theRapidity) > JpsiYMin && fabs(theRapidity) < JpsiYMax &&
	ok2 && ok1 &&
	(trig == 1 || trig2 == 1) &&
	vprob >0.01 &&
	Jq == 0){
      
      //cout << Jq << endl;
      Jpsi_Pt->setVal(thePt); 
      Jpsi_Y->setVal(theRapidity); 
      Jpsi_Mass->setVal(theMass);
      Psip_Mass->setVal(theMass);
      Jpsi_Ct->setVal(theCt);
      Jpsi_CtErr->setVal(theCtErr);
      Jpsi_Type->setIndex(theCat,kTRUE);
      //cout << "Category " << theCat << endl;

      RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr);
      RooArgList varlist2_tmp(*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr);
      RooArgList varlist3_tmp(*Jpsi_Mass,*Psip_Mass,*Jpsi_Pt,*Jpsi_Y,*Jpsi_Type,*Jpsi_Ct,*Jpsi_CtErr,*Jpsi_PsiP);

      for (unsigned int i = 0; i < rapRegions; i++) {
	if (fabs(theRapidity) < rapLimits[i+1] && fabs(theRapidity) > rapLimits[i]) {
	  if (doMerge) {
	    if (doWideRange) {
	      dataJpsi[i]->add(varlist_tmp);
	    } else {
	      if (theMass <= 3.3) {
		Jpsi_PsiP->setIndex(0,kTRUE); dataJpsi[i]->add(varlist3_tmp);
	      } else if (theMass > 3.3 && theMass < 3.5) {
		Jpsi_PsiP->setIndex(0,kTRUE); dataJpsi[i]->add(varlist3_tmp);
		Jpsi_PsiP->setIndex(1,kTRUE); dataJpsi[i]->add(varlist3_tmp);
	      } else {
		Jpsi_PsiP->setIndex(1,kTRUE); dataJpsi[i]->add(varlist3_tmp);
	      }
	    }
	  } else {
	    if (theMass < 3.5) dataJpsi[i]->add(varlist_tmp);
	    if (theMass > 3.3) dataPsip[i]->add(varlist2_tmp);
	  }
	}
      }
    }
  }

  TFile* Out[rapRegions];
  char namefile[200];
  for (unsigned int i = 0; i < rapRegions; i++) {
    if (doMerge && !doWideRange) {
      sprintf(namefile,"datasets/DataMerged2010_rap%d-%d.root",int(rapLimits[i]*10),int(rapLimits[i+1]*10));
    } else if (doMerge && doWideRange) {
      sprintf(namefile,"datasets/DataMergedWide2010_rap%d-%d.root",int(rapLimits[i]*10),int(rapLimits[i+1]*10)); 
    } else {
      sprintf(namefile,"datasets/Data2010_rap%d-%d.root",int(rapLimits[i]*10),int(rapLimits[i+1]*10));
    }
    Out[i] = new TFile(namefile,"RECREATE");
    Out[i]->cd();
    dataJpsi[i]->Write();
    if (!doMerge) dataPsip[i]->Write();
    Out[i]->Close();
  }

  // delete dataJpsi;
  // delete dataPsip;
  delete JP;
  delete m1P;
  delete m2P;
  delete Jpsi_Mass;
  delete Jpsi_Pt;
  delete Jpsi_Y;
  delete file;
}


