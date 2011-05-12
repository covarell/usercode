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
#include <TLatex.h>

#include "RooFit.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include <TCanvas.h>
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TVector3.h>
#include <math.h>

using namespace RooFit;
using namespace std;

void getrange(string &varRange, float *varmin, float *varmax)
{
 if (sscanf(varRange.c_str(), "%f-%f", varmin, varmax) == 0) {
   cout << varRange.c_str() << ": range not valid!" << endl;
   assert(0);
  }

 return;
}

bool Cowboy(int mu1_charge, TLorentzVector* mu1, TLorentzVector *mu2){
  return (mu1_charge * mu1->DeltaPhi(*mu2) >0.);
} 

bool isAccept(const TLorentzVector* aMuon) {
   // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
   return (fabs(aMuon->Eta()) < 2.4 &&
           ((fabs(aMuon->Eta()) < 1.3 && aMuon->Pt() > 3.3) ||
           (fabs(aMuon->Eta()) > 1.3 && fabs(aMuon->Eta()) < 2.2 && aMuon->P() > 2.9) ||
           (fabs(aMuon->Eta()) > 2.2 && aMuon->Pt() > 0.8)));

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

void defineMassSignal(RooWorkspace *ws)
{
  //SIGNAL FUNCTION CANDIDATES:

  //Normal Gaussians
  ws->factory("Gaussian::signalG1(Jpsi_Mass,meanSig1[3.0975,3.05,3.15],sigmaSig1[0.02,0.008,0.2])");
  ws->factory("Gaussian::signalG2(Jpsi_Mass,meanSig2[3.0975,3.05,3.15],sigmaSig2[0.03,0.008,0.2])");

  //Gaussian with same mean as signalG1
  ws->factory("Gaussian::signalG2OneMean(Jpsi_Mass,meanSig1,sigmaSig2)");

  //Crystall Ball
  ws->factory("CBShape::sigCB(Jpsi_Mass,meanSig1,sigmaSig2,alpha[0.5,0.,3.],enne[5.,1.,30.])");

  //SUM OF SIGNAL FUNCTIONS

  //Sum of Gaussians with different mean
  // ws->factory("SUM::sigPDF(coeffGauss[0.1,0.,1.]*signalG1,signalG2)");

  //Sum of Gaussians with same mean
  ws->factory("SUM::sigPDFOneMean(coeffGauss[0.1,0.,1.]*signalG1,signalG2OneMean)");

  //Sum of a Gaussian and a CrystalBall
  // ws->factory("SUM::sigCBGauss(coeffGauss*sigCB,signalG2)");

  //Sum of a Gaussian and a CrystalBall
  ws->factory("SUM::sigCBGaussOneMean(coeffGauss*sigCB,signalG1)");

  return;
}

void defineMassBackground(RooWorkspace *ws)
{
  //Second order polynomial, the 2nd coefficient is by default set to zero
  ws->factory("Polynomial::CPolFunct(Jpsi_Mass,{CoefPol1[-0.05,-1500.,1500.],CoefPol2[-1.,-10.,0.]})");

  //Exponentials
  ws->factory("Exponential::expFunctP(Jpsi_Mass,coefExpP[-1.,-3.,-0.6])");
  ws->factory("Exponential::expFunctF(Jpsi_Mass,coefExpF[-1.,-3.,-0.6])");
  return;
}

int main(int argc, char* argv[]) {
  
  float JpsiMassMin=2.6;
  float JpsiMassMax=3.5;
  float JpsiPtMin=0;
  float JpsiPtMax=50;
  float JpsiYMin=0;
  float JpsiYMax=2.4;
  float JpsiCtMin = -2.0;
  float JpsiCtMax = 3.5;

  char fileName[100];
  string prange, yrange;	
  unsigned int isTree;  

  if ( argc < 5 ){
    cout << "missing arguments..." << endl; 
    return 1;
  }
  strcpy(fileName,argv[1]);
  prange = argv[2];
  yrange = argv[3];
  isTree = atoi(argv[4]);

  getrange(prange,&JpsiPtMin,&JpsiPtMax);
  getrange(yrange,&JpsiYMin,&JpsiYMax);

  cout << "pT = " << prange << ", y = " << yrange << endl;
  cout << "pT = " << JpsiPtMin << "-" << JpsiPtMax << ", y = " << JpsiYMin << "-" << JpsiYMax << endl; 

  TFile *file= TFile::Open(fileName);
  RooWorkspace *ws = new RooWorkspace("ws");  

  TLorentzVector* JP= new TLorentzVector;
  TLorentzVector* m1P= new TLorentzVector;
  TLorentzVector* m2P= new TLorentzVector;
  
  float vprob, theCt, theCtErr;
  int trig,trig2,theCat,Jq;

  RooDataSet* dataJpsi;
  RooRealVar* Jpsi_Mass;
  RooRealVar* Jpsi_Pt;
  RooRealVar* Jpsi_Ct;
  RooRealVar* Jpsi_CtErr;
  RooRealVar* Jpsi_Y;
  RooCategory* Jpsi_Type;
  RooCategory* Pass_Fail;

  Jpsi_Mass = new RooRealVar("Jpsi_Mass","J/psi mass",JpsiMassMin,3.5,"GeV/c^{2}");
  Jpsi_Pt = new RooRealVar("Jpsi_Pt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
  Jpsi_Y = new RooRealVar("Jpsi_Y","J/psi y",-JpsiYMax,JpsiYMax);
  Jpsi_Type = new RooCategory("Jpsi_Type","Category of Jpsi");
  Jpsi_Ct = new RooRealVar("Jpsi_Ct","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
  Jpsi_CtErr = new RooRealVar("Jpsi_CtErr","J/psi ctau error",0.,4.,"mm");
  Pass_Fail = new RooCategory("Pass_Fail","Pass or fail");
  
  Jpsi_Mass->setBins(20);
  ws->import(*Jpsi_Mass);
  
  Jpsi_Type->defineType("GG",0);
  Jpsi_Type->defineType("GT",1);
  Jpsi_Type->defineType("TT",2);

  Pass_Fail->defineType("fail",0);
  Pass_Fail->defineType("pass",1);

  RooArgList varlist(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Pass_Fail);

  if (isTree) {

    TTree * Tree=(TTree*)file->Get("data");

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

    int n = Tree->GetEntries();
    int nev=n;
    dataJpsi = new RooDataSet("dataJpsi","A sample",varlist);

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

      //cout << "Accepted1 " << ok1 << "  Accepted2 " << ok2 << endl;
      if (theMass > JpsiMassMin && theMass < JpsiMassMax && 
	theCt > JpsiCtMin && theCt < JpsiCtMax && 
	thePt > JpsiPtMin && thePt < JpsiPtMax && 
	fabs(theRapidity) > JpsiYMin && fabs(theRapidity) < JpsiYMax &&
	ok2 && ok1 &&
	vprob >0.001 &&
	Jq == 0){
      
        //cout << Jq << endl;
        Jpsi_Pt->setVal(thePt); 
        Jpsi_Y->setVal(theRapidity); 
        Jpsi_Mass->setVal(theMass);
        Jpsi_Ct->setVal(theCt);
        Jpsi_CtErr->setVal(theCtErr);
        Jpsi_Type->setIndex(theCat,kTRUE);
        //cout << "Category " << theCat << endl;

        if (trig == 1 || trig2 == 1) {
	  Pass_Fail->setIndex(1,kTRUE);
        } else {
          Pass_Fail->setIndex(0,kTRUE);
        }

        RooArgList varlist_tmp(*Jpsi_Mass,*Jpsi_Pt,*Jpsi_Y,*Pass_Fail);
        dataJpsi->add(varlist_tmp);
      }
    }

    TFile* Out;
    char namefile[200];
    sprintf(namefile,"datasets/multijet_pT%s_y%s.root",prange.c_str(),yrange.c_str());
    Out = new TFile(namefile,"RECREATE");
    Out->cd();
    dataJpsi->Write();
    Out->Close();
  
  } else {

    string titlestr;
    char theCut[200];
    sprintf(theCut,"Jpsi_Pt < %f && Jpsi_Pt > %f && Jpsi_Y < %f && Jpsi_Y > %f",JpsiPtMax,JpsiPtMin,JpsiYMax,JpsiYMin); 
    dataJpsi = (RooDataSet*)file->Get("dataJpsi");
    // FIT
    defineMassSignal(ws);
    defineMassBackground(ws);

    ws->var("alpha")->setVal(1.523);  
    ws->var("coefExpP")->setVal(-1.166);
    ws->var("coefExpF")->setVal(-1.166); 
    ws->var("coeffGauss")->setVal(0.4664);
    ws->var("enne")->setVal(1.394);   
    ws->var("meanSig1")->setVal(3.096); 
    ws->var("sigmaSig1")->setVal(0.0311); 
    ws->var("sigmaSig2")->setVal(0.01895);  

    ws->var("alpha")->setConstant(kTRUE);
    ws->var("coeffGauss")->setConstant(kTRUE);
    ws->var("enne")->setConstant(kTRUE);
    ws->var("meanSig1")->setConstant(kTRUE);
    // ws->var("sigmaSig1")->setConstant(kTRUE);
    ws->var("sigmaSig2")->setConstant(kTRUE);

    RooDataSet *reddata = (RooDataSet*) dataJpsi->reduce(theCut);
    RooDataSet *reddataP = (RooDataSet*) reddata->reduce("Pass_Fail == Pass_Fail::pass");
    RooDataSet *reddataF = (RooDataSet*) reddata->reduce("Pass_Fail == Pass_Fail::fail");

    ws->factory("SUM::massPDFP(NSigP[5000.,2.,10000000.]*sigCBGaussOneMean,NBkgP[2000.,2.,10000000.]*expFunctP)");
    RooRealVar eff("eff","#epsilon",0.75,0.001,1.);
    ws->import(eff);
    RooFormulaVar NSigF("NSigF", "@0*(1./@1 - 1.)", RooArgList(*(ws->var("NSigP")), *(ws->var("eff"))));  
    ws->import(NSigF);  
    ws->factory("SUM::massPDFF(NSigF*sigCBGaussOneMean,NBkgF[2000.,2.,10000000.]*expFunctF)");
    // ws->pdf("massPDFP")->fitTo(*reddataP,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2));
    // ws->pdf("massPDFF")->fitTo(*reddataF,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2));
    RooSimultaneous massSim("massSim","mass simultaneous PDF",RooArgList(*(ws->pdf("massPDFF")),*(ws->pdf("massPDFP"))),*Pass_Fail);
    ws->import(massSim);
    ws->pdf("massSim")->fitTo(*reddata,Extended(1),Minos(0),SumW2Error(kTRUE),NumCPU(2)); 

    gROOT->ProcessLine(".L mytdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    ws->var("Jpsi_Mass")->SetTitle("J/#psi mass");

    RooPlot *mframe = ws->var("Jpsi_Mass")->frame();

    titlestr = "Passing dimuons, p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    mframe->SetTitle(titlestr.c_str());

    reddata->plotOn(mframe,DataError(RooAbsData::SumW2),RooFit::Cut("Pass_Fail == Pass_Fail::pass"));

    ws->pdf("massPDFP")->plotOn(mframe,LineColor(kRed),Components("expFunctP"),LineStyle(kDashed),Normalization(reddataP->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("massPDFP")->plotOn(mframe,LineColor(kBlue),Normalization(reddataP->sumEntries(),RooAbsReal::NumEvent));
    
    RooPlot *mframe2 = ws->var("Jpsi_Mass")->frame();

    titlestr = "Failing dimuons, p_{T} = " + prange + " GeV/c and |y| = " + yrange;
    mframe2->SetTitle(titlestr.c_str());

    reddata->plotOn(mframe2,DataError(RooAbsData::SumW2),RooFit::Cut("Pass_Fail == Pass_Fail::fail"));

    ws->pdf("massPDFF")->plotOn(mframe2,LineColor(kRed),Components("expFunctF"),LineStyle(kDashed),Normalization(reddataF->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("massPDFF")->plotOn(mframe2,LineColor(kBlue),Normalization(reddataF->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("massPDFF")->paramOn(mframe2,Parameters(RooArgSet(eff)));

    TCanvas c1("c1","c1",10,10,800,400);
    c1.Divide(2,1);
    c1.cd(1);   mframe->Draw();
    c1.cd(2);   mframe2->Draw();

    TLatex *t = new TLatex();
    t->SetNDC();
    t->SetTextAlign(22);
    t->SetTextSize(0.035);
    titlestr = "pT = " + prange;
    t->DrawLatex(0.8,0.84,titlestr.c_str()); 
    titlestr = "|y| = " + yrange;
    t->DrawLatex(0.8,0.78,titlestr.c_str()); 
    c1.SaveAs("/afs/cern.ch/user/c/covarell/public/html/quarkonia/testEff.gif");
    
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


