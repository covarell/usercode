#include "RooRealVar.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <sstream>
#include <vector>
#include "../src/AngularPdfFactory.cc"
#include "../PDFs/RooqqZZ_JHU.h"
#include "../PDFs/RooTsallis.h"

/* - - - - - - - - - - - - - - - - - - - - - - - 
================================================

be sure to compile/load everything:

root -l -n 
gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L MELA.C+
 - -  - - - - - - - - - - - - - - -  - - -  -  -
===============================================*/

using namespace RooFit ;

const int mZZbins=350;
int lowMzz=100;
int highMzz=800;
int lowM2=12;

double binning[mZZbins] = {
  100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 
  120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 
  140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 
  160, 162, 164, 166, 168, 170, 172, 174, 176, 178,
  180, 182, 184, 186, 188, 190, 192, 194, 196, 198,
  200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 
  220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 
  240, 242, 244, 246, 248, 250, 255, 260, 265, 270, 
  275, 280, 285, 290, 295, 300, 310, 320, 330, 340, 
  350, 360, 370, 380, 390, 400, 420, 440, 460, 480, 
  500, 520, 560, 600, 640, 680, 720, 760, 800 };

// use for variable binning
//TH1F* mzzBinning = new TH1F("mzzBinning","mzzBinning",108,binning);
// use for constant binning
TH1F* mzzBinning = new TH1F("mzzBinning","mzzBinning",350,100,800);

TFile *tempf = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");


template <typename U>
void checkZorder(U& z1mass, U& z2mass,
                 U& costhetastar, U& costheta1,
                 U& costheta2, U& phi,
                 U& phistar1, U& pt, U& ipsilon){

  U tempZ1mass=z1mass;
  U tempZ2mass=z2mass;
  U tempH1=costheta1;
  U tempH2=costheta2;
  U tempHs=costhetastar;
  U tempPhi1=phistar1;
  U tempPhi=phi;
  U tempPt=pt;
  U tempY=ipsilon;

  if(z2mass>z1mass){

    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    pt=tempPt;
    ipsilon=tempY;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*3.1415;
    if(phistar1<-3.1415)
      phistar1=phistar1+2*3.1415;

  }else
    return;

}

vector<double> my8DTemplate(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1){

  //read from a file the 3D and 2D template
  TH1F *h_mzz= (TH1F*)(tempf->Get("h_mzz"));
  TH3F *h_mzzm1m2= (TH3F*)(tempf->Get("h_mzzm1m2"));
  TH2F *h_mzzcosthetastar= (TH2F*)(tempf->Get("h_mzzcosthetastar"));
  TH2F *h_mzzcostheta1= (TH2F*)(tempf->Get("h_mzzcostheta1"));
  TH2F *h_mzzcostheta2= (TH2F*)(tempf->Get("h_mzzcostheta2"));
  TH2F *h_mzzphi1= (TH2F*)(tempf->Get("h_mzzphi1"));
  TH2F *h_mzzphi= (TH2F*)(tempf->Get("h_mzzphi"));

  //multiply the P values
  double n = h_mzz->GetBinContent(h_mzz->FindBin(mZZ));
  double Pmzzm1m2 = h_mzzm1m2->GetBinContent(h_mzzm1m2->FindBin(mZZ,m1,m2));

  // - - - - - - - - - - - - - - - whitbeck
  // if bin has no events: add 1
  // safety feature to prevent LD = 1 as a
  // result of low statistics

  if(Pmzzm1m2==0){
    Pmzzm1m2++;
    }
  // - - - - - - - - - - - - - - - 

  double Pmzzcosthetastar = h_mzzcosthetastar->GetBinContent(h_mzzcosthetastar->FindBin(mZZ,costhetastar));
  double Pmzzcostheta2 = h_mzzcostheta2->GetBinContent(h_mzzcostheta2->FindBin(mZZ,costheta2));
  double Pmzzcostheta1 = h_mzzcostheta1->GetBinContent(h_mzzcostheta1->FindBin(mZZ,costheta1));
  double Pmzzphi1 = h_mzzphi1->GetBinContent(h_mzzphi1->FindBin(mZZ,phi1));
  double Pmzzphi = h_mzzphi->GetBinContent(h_mzzphi->FindBin(mZZ,phi));

  //normalization
  double binwidth_mzzm1m2 = h_mzzm1m2->GetYaxis()->GetBinWidth(1) * h_mzzm1m2->GetZaxis()->GetBinWidth(1);
  double binwidth_mzzcosthetastar = h_mzzcosthetastar->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzcostheta1 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzcostheta2 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzphi1 = h_mzzphi1->GetYaxis()->GetBinWidth(1);
  double binwidth_mzzphi = h_mzzphi->GetYaxis()->GetBinWidth(1);

  double Pmzzm1m2_norm = Pmzzm1m2/(n*binwidth_mzzm1m2); 
  double Pmzzcosthetastar_norm = Pmzzcosthetastar/(n*binwidth_mzzcosthetastar);
  double Pmzzcostheta1_norm = Pmzzcostheta1/(n*binwidth_mzzcostheta1);
  double Pmzzcostheta2_norm = Pmzzcostheta2/(n*binwidth_mzzcostheta2);
  double Pmzzphi1_norm = Pmzzphi1/(n*binwidth_mzzphi1);
  double Pmzzphi_norm = Pmzzphi/(n*binwidth_mzzphi);

  vector <double> P;
  P.push_back(Pmzzm1m2);
  P.push_back(Pmzzcosthetastar);
  P.push_back(Pmzzcostheta1);
  P.push_back(Pmzzcostheta2);
  P.push_back(Pmzzphi);
  P.push_back(Pmzzphi1);

  vector <double> P_norm;
  P_norm.push_back(Pmzzm1m2_norm);
  P_norm.push_back(Pmzzcosthetastar_norm);
  P_norm.push_back(Pmzzcostheta1_norm);
  P_norm.push_back(Pmzzcostheta2_norm);
  P_norm.push_back(Pmzzphi_norm);
  P_norm.push_back(Pmzzphi1_norm);

  if(normalized)
    return P_norm;
  else
    return P;
}

//=======================================================================

pair<double,double> likelihoodDiscriminant (double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1,bool withPtY = false, double pt = 0.0, double absy = 0.0, double scaleFactor=5.0){

  RooRealVar* z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0,180);
  RooRealVar* z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0,120); 
  RooRealVar* costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1,1);  
  RooRealVar* costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
  RooRealVar* phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  RooRealVar* costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1,1); 
  RooRealVar* phi1_rrv= new RooRealVar("phi1","#Phi^{*}_{1}",-3.1415,3.1415);
  RooRealVar* mzz_rrv= new RooRealVar("mzz","mZZ",80,1000);
  RooRealVar* absy_rrv= new RooRealVar("absy","|y_{ZZ}|",0.0,4.0);
  RooRealVar* pt_rrv= new RooRealVar("pt","p_{T,ZZ}",0.0,1000.);

  static const int Nptparams = 11;

  char* rrvnamesB[Nptparams] = {"m","n0","n1","n2","ndue","bb0","bb1","bb2","T0","T1","T2"};
  RooRealVar *ptparamsB[Nptparams];
  RooArgSet* allparamsB = new RooArgSet();
  for (int i = 0; i < Nptparams; i++) {
    ptparamsB[i] = new RooRealVar(rrvnamesB[i],rrvnamesB[i],-10000.,10000.);
    allparamsB->add(*ptparamsB[i]);
  }

  char* rrvnamesS[Nptparams] = {"ms","ns0","ns1","ns2","ndues","bbs0","bbs1","bbs2","Ts0","Ts1","Ts2"};
  RooRealVar *ptparamsS[Nptparams];
  RooArgSet* allparamsS = new RooArgSet();
  for (int i = 0; i < Nptparams; i++) {
    ptparamsS[i] = new RooRealVar(rrvnamesS[i],rrvnamesS[i],-10000.,10000.);
    allparamsS->add(*ptparamsS[i]);
  }

  AngularPdfFactory *SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);
  RooqqZZ_JHU* SMZZ = new RooqqZZ_JHU("SMZZ","SMZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv);

  RooTsallis* sigPt = new RooTsallis("sigPt","sigPt",*pt_rrv,*mzz_rrv,
				     *ptparamsS[0],*ptparamsS[1],*ptparamsS[2],
				     *ptparamsS[3],*ptparamsS[4],*ptparamsS[5],
				     *ptparamsS[6],*ptparamsS[7],*ptparamsS[8],
				     *ptparamsS[9],*ptparamsS[10]);
  if (withPtY) allparamsS->readFromFile("../datafiles/allParamsSig.txt",0);

  RooTsallis* bkgPt = new RooTsallis("bkgPt","bkgPt",*pt_rrv,*mzz_rrv,
				     *ptparamsB[0],*ptparamsB[1],*ptparamsB[2],
				     *ptparamsB[3],*ptparamsB[4],*ptparamsB[5],
				     *ptparamsB[6],*ptparamsB[7],*ptparamsB[8],
				     *ptparamsB[9],*ptparamsB[10]);
  if (withPtY) allparamsB->readFromFile("../datafiles/allParamsBkg.txt",0);

  checkZorder<double>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1,pt,absy);

  z1mass_rrv->setVal(m1);  
  z2mass_rrv->setVal(m2);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  costhetastar_rrv->setVal(costhetastar);
  phi1_rrv->setVal(phi1);
  mzz_rrv->setVal(mZZ);
  pt_rrv->setVal(pt);
  absy_rrv->setVal(absy);

  vector <double> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phi1);

  double Pbackg=0.;
  double Psig=0.;

  if(mZZ>80 && mZZ<180){
    Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*5.0;
    Psig = SMHiggs->getVal(mZZ);
  }if(mZZ>180&&mZZ<=2*91.188){
    z1mass_rrv->setVal(mZZ/2. - 1e-9);  // turns out signal norm will by zero if m1+m2==mzz
    z2mass_rrv->setVal(mZZ/2. - 1e-9);  // adding tiny amount to avoid this
    Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
  }if(mZZ>2*91.188){
    z1mass_rrv->setVal(91.188);
    z2mass_rrv->setVal(91.188);
    Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
    Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
  }

  if (withPtY) {
    Pbackg *= bkgPt->getVal()/(bkgPt->createIntegral(RooArgSet(*pt_rrv))->getVal());
    Psig *= sigPt->getVal()/(sigPt->createIntegral(RooArgSet(*pt_rrv))->getVal());
  }

  // - - - - - - - - - - - - - - - - - - - - - Whitbeck 
  // check whether P[i] is zero and print warning
  // message if so

  char* varName[6]={"m1/m2","costhetastar","costheta1","coshteta2","phi","phi1"};
  for(int iVar=0; iVar<6; iVar++){

    if(P[iVar]==0 && (m1+m2)<mZZ && m2>4 && mZZ>80 && mZZ<180)
	cout << " uh oh... Probability of " << varName[iVar] << " is zero." << endl;
  }
  // - - - - - - - - - - - - - - - - - - - - - 
 
  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete costhetastar_rrv;
  delete phi1_rrv;
  delete mzz_rrv;
  delete pt_rrv;
  delete absy_rrv;
  delete sigPt;
  delete bkgPt;
  delete SMZZ;
  delete SMHiggs;

  return make_pair(Psig,Pbackg);
}

//=======================================================================

void addDtoTree(char* inputFile, float minMzz = 100., float maxMzz = 1000., bool containsPtY = false){

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

  char inputFileName[100];
  char outputFileName[150];
  sprintf(inputFileName,"%s.root",inputFile);

  sprintf(outputFileName,"%s_%d-%d_withDiscriminants.root",inputFile,int(minMzz),int(maxMzz));

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
    if(sigFile)
        sigTree = (TTree*) sigFile->Get("SelectedTree");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }
  
 
  // sprintf(outputFileName,"signal_withDiscriminants.root",inputFile);

  // sprintf(outputFileName,"background_withDiscriminants.root",inputFile);

  // TChain* sigTree = new TChain("SelectedTree");
  /* sigTree->Add("../datafiles/4e/HZZ4lTree_H125.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_H125.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_H125.root"); */
  // sigTree->Add("../datafiles/4e/HZZ4lTree_ggZZ2l2l.root");
  // sigTree->Add("../datafiles/4e/HZZ4lTree_ggZZ4l.root");
  /* sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo2e2mu.root");
  sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo2e2tau.root");
  sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo2mu2tau.root");
  sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4e.root");
  sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4mu.root");
  sigTree->Add("../datafiles/4e/HZZ4lTree_ZZTo4tau.root");
  // sigTree->Add("../datafiles/4mu/HZZ4lTree_ggZZ2l2l.root");
  // sigTree->Add("../datafiles/4mu/HZZ4lTree_ggZZ4l.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo2e2mu.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo2e2tau.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo2mu2tau.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4e.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4mu.root");
  sigTree->Add("../datafiles/4mu/HZZ4lTree_ZZTo4tau.root");
  // sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ggZZ2l2l.root");
  // sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ggZZ4l.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2mu.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2e2tau.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo2mu2tau.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo4e.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo4mu.root");
  sigTree->Add("../datafiles/2mu2e/HZZ4lTree_ZZTo4tau.root"); */

  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = new TTree("newTree","SelectedTree"); 

  float m1,m2,mzz,h1,h2,hs,phi,phi1,D;
  float ZZPt, ZZY, Dpt, w;
  
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  // sigTree->SetBranchAddress("costheta1",&h1); 
  // sigTree->SetBranchAddress("costheta2",&h2);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1); 
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  // sigTree->SetBranchAddress("phi",&phi);  
  // sigTree->SetBranchAddress("phistar1",&phi1);
  sigTree->SetBranchAddress("helphi",&phi);  
  sigTree->SetBranchAddress("phistarZ1",&phi1);

  if (containsPtY) {
    sigTree->SetBranchAddress("ZZPt",&ZZPt);
    // newTree->SetBranchAddress("ZZY",&ZZY);
  }
  sigTree->SetBranchAddress("MC_weight",&w);

  newTree->Branch("z1mass",&m1,"z1mass/F");
  newTree->Branch("z2mass",&m2,"z2mass/F");
  newTree->Branch("zzmass",&mzz,"zzmass/F");
  newTree->Branch("costheta1",&h1,"costheta1/F"); 
  newTree->Branch("costheta2",&h2,"costheta2/F");
  newTree->Branch("costhetastar",&hs,"costhetastar/F");
  newTree->Branch("phi",&phi,"phi/F");  
  newTree->Branch("phistar1",&phi1,"phistar1/F");
  if (containsPtY) {
    newTree->Branch("ZZPt",&ZZPt,"ZZPt/F");
    // newTree->Branch("ZZY",&ZZY,"ZZY/F");
    newTree->Branch("melaLDWithPt",&Dpt,"melaLDWithPt/F");
  }
  newTree->Branch("melaLD",&D,"melaLD/F");
  newTree->Branch("MC_weight",&w,"MC_weight/F");

  for(int iEvt=0; iEvt<sigTree->GetEntries(); iEvt++){

    if(iEvt%5000==0) 
      cout << "event: " << iEvt << endl;

    sigTree->GetEntry(iEvt);

    checkZorder<float>(m1,m2,hs,h1,h2,phi,phi1,ZZPt,ZZY);

    if(mzz>minMzz && mzz<maxMzz && m2>12. ) 
      {

      //MELA LD
	pair<double,double> P;
	P = likelihoodDiscriminant(mzz, m1, m2, hs, h1, h2, phi, phi1);
	D=P.first/(P.first+P.second);

	if (containsPtY) {
	  pair<double,double> P2;
	  P2 = likelihoodDiscriminant(mzz, m1, m2, hs, h1, h2, phi, phi1, true, ZZPt, fabs(ZZY));
 
	  Dpt=P2.first/(P2.first+P2.second);
	}

      newTree->Fill();

    }
   }

  newFile->cd();
  newTree->Write("angles"); 
  newFile->Close();

}

//=======================================================================

vector<TH1F*> LDDistributionBackground(TTree* chain){

  float mZZ, m2, m1, costhetastar, costheta1, costheta2, phi, phi1, pt, y;
  float MC_weight=1;
  float mela=-99;
  chain->SetBranchAddress("ZZMass",&mZZ);
  chain->SetBranchAddress("MC_weight_noxsec",&MC_weight);  
  chain->SetBranchAddress("ZZLD",&mela);

  //TFile *f = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

  TH1F *h_LDbackground= new TH1F("LD_background","LD_background",31,0,1.01);
  vector<TH1F*> vh_LDbackground;

  // use for variable binning
  //for (int i=1; i<mZZbins; i++){
  // use for constant binning 
  for (int i=1; i<=mZZbins; i++){
    std::string s;
    std::stringstream out;
    out << i;
    s = out.str();
    TString name = "h_LDbackground_"+s;
    vh_LDbackground.push_back((new TH1F(name,name,30,0,1)));
    vh_LDbackground[i-1]->Sumw2();
  }

  for (Int_t i=0; i<chain->GetEntries();i++) {

    mela=-99;

    if(i%100000 ==0)
      cout<<"event "<<i<<endl;
      
    chain->GetEvent(i); 

    if( !(mZZ>lowMzz && mZZ<highMzz) ) // && m2>lowM2) )
      continue;
          
    if(mela<0){

      checkZorder<float>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1,pt,y);
      
      pair<double,double> P =  likelihoodDiscriminant(mZZ, m1, m2, costhetastar, costheta1, costheta2, phi, phi1);
      
      mela=P.first/(P.first+P.second);
    }
      
    h_LDbackground->Fill(mela,MC_weight);
    (vh_LDbackground[mzzBinning->FindBin(mZZ)-1])->Fill(mela,MC_weight);

  }

  TCanvas *LD = new TCanvas("LD_background","LD_background",400,400);
  h_LDbackground->Draw();
  LD->Print("LD_background.eps");

  vh_LDbackground.push_back(h_LDbackground);
 
  return vh_LDbackground;
}


//=======================================================================
vector<TH1F*> LDDistributionSignal(TTree* chain){
 
  float mZZ, m2, m1, costhetastar, costheta1, costheta2, phi, phi1, pt, y;
  float MC_weight=1;
  float mela=-99;
  chain->SetBranchAddress("ZZMass",&mZZ);
  chain->SetBranchAddress("MC_weight_noxsec",&MC_weight);
  chain->SetBranchAddress("ZZLD",&mela);

  //TFile *f = new TFile("../datafiles/my8DTemplateNotNorm.root","READ");

  TH1F *h_LDsignal= new TH1F("LD_signal","LD_signal",31,0,1.01);
  vector<TH1F*> vh_LDsignal;

  // use for variable binning
  //for (int i=1; i<mZZbins; i++){
  // use for constant binning
  for (int i=1; i<=mZZbins; i++){
    std::string s;
     std::stringstream out;
     out << i;
     s = out.str();
     TString name = "  %h_LDsignal_"+s;
     vh_LDsignal.push_back((new TH1F(name,name,30,0,1)));
     vh_LDsignal[i-1]->Sumw2();
  }

  for (Int_t i=0; i<chain->GetEntries();i++) {

    mela=-99;

    if(i%100000 ==0)
      cout<<"event "<<i<<endl;
    
    chain->GetEvent(i); 

    if( !(mZZ>lowMzz && mZZ<highMzz)) // && m2>lowM2) )
      continue;
      
    if(mela<0){

      checkZorder<float>(m1,m2,costhetastar,costheta1,costheta2,phi,phi1,pt,y);
      
      pair<double,double> P =  likelihoodDiscriminant(mZZ, m1, m2, costhetastar, costheta1, costheta2, phi, phi1);

      mela=P.first/(P.first+P.second);

    }

    h_LDsignal->Fill(mela,MC_weight);
    (vh_LDsignal[mzzBinning->FindBin(mZZ)-1])->Fill(mela,MC_weight);
     
  }
    
  cout << "Drawing" << endl;

  TCanvas *LD = new TCanvas("LD_signal","LD_signal",400,400);
  h_LDsignal->Draw();
  LD->Print("LD_signal.eps");

  cout << "push_back" << endl;

  vh_LDsignal.push_back(h_LDsignal);

  cout << "done" << endl;

  return vh_LDsignal;
}

//=======================================================================

TH2F* smoothTemplate(TH2F* oldTemp, TH2F* numEvents){

  int effectiveArea=3;

  TH2F* newTemp = new TH2F(*oldTemp);

  double average=0;
  int nBins=0;

  for(int i=1; i<mZZbins; i++){
    
    for(int j=1; j<=30; j++){
      
      if( numEvents->GetBinContent(i,j)!=0 && numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)<.5 )  continue;

      if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i>40 && i<100 && j>15){
	effectiveArea=2;
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}

	if(average==0){
	  effectiveArea=3;
	  nBins=0;

	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }
	}
	if(average==0){
	  effectiveArea=4;
	  nBins=0;

	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  } 
	}
      }else if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i>100 ){
	effectiveArea=2;

	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}

	if(average==0){
	  effectiveArea=3;
	  nBins=0;
	  
	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }
	}
	if(average==0){
	  effectiveArea=4;
	  nBins=0;
	  
	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<41 || i+a>mZZbins-1 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }

	}

      }else if( (numEvents->GetBinContent(i,j)==0 || numEvents->GetBinError(i,j)/numEvents->GetBinContent(i,j)>.5)&& i<=40 ){
	effectiveArea=1;
	
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  for(int b=-effectiveArea; b<=effectiveArea; b++){
	    if( i+a<1 || i+a>40 || j+b<1 || j+b>30 ) continue;
	    average+=oldTemp->GetBinContent(i+a,j+b);
	    nBins++;
	  }
	}

	if(average==0){
	  effectiveArea=2;
	  nBins=0;

	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<1 || i+a>40 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }
	}

	if(average==0){
	  effectiveArea=3;
	  nBins=0;

	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<1 || i+a>40 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }
	}

	if(average==0){
	  effectiveArea=4;
	  nBins=0;

	  for(int a=-effectiveArea; a<=effectiveArea; a++){
	    for(int b=-effectiveArea; b<=effectiveArea; b++){
	      if( i+a<1 || i+a>40 || j+b<1 || j+b>30 ) continue;
	      average+=oldTemp->GetBinContent(i+a,j+b);
	      nBins++;
	    }
	  }
	}

      }else continue;

      newTemp->SetBinContent(i,j,average/nBins);
      average=0;
      nBins=0;
    }// end loop over D bins
    //if(i<40) cout << endl;
  }// end loop over mZZ bins

  double norm=0;

  for(int i=1; i<=mZZbins; i++){
    for(int j=1; j<=30; j++){
      norm+=newTemp->GetBinContent(i,j);
    }

    for(int j=1; j<=30; j++){
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);
    }

    norm=0;

  }

  return newTemp;

}

//=======================================================================
pair<TH2F*,TH2F*> reweightForCRunc(TH2F* temp){

  cout << "reweightForCRunc" << endl;

  TH2F* tempUp = new TH2F(*temp);
  TH2F* tempDn = new TH2F(*temp);

  pair<TH2F*,TH2F*> histoPair(0,0);

  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;
  int point=-1;

  const int numPoints=8;

  double low[numPoints]   ={100.,        120.,        140.,         160.,     180.,     220.,     260.,     300. }; 
  double high[numPoints]  ={120.,        140.,        160.,         180.,     220.,     260.,     300.,     1000.};
  double slope[numPoints] ={4.71836e-01, 1.17671e-01, -3.81680e-01, -1.20481, -1.21944, -2.06928, -1.35337, 0.0 };
  double yIntr[numPoints] ={6.83860e-01, 9.38454e-01, 1.12690,      1.24502,  1.72764,  2.11050,  1.52771,  1.0 }; 

  for(int i=1; i<=temp->GetNbinsX(); i++){
    point = -1;

    // choose correct scale factor
    for(int p=0; p<numPoints; p++){
      if( (i*2.+101.)>=low[p] && (i*2.+101.)<high[p] ){
	point = p;
      }
    }
    if(point == -1){
      cout << "ERROR: could not find correct scale factor"<< endl;
      return histoPair;
    }

    for(int j=1; j<=temp->GetNbinsY(); j++){

      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(slope[point]*(double)j/30.+yIntr[point]);
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(-slope[point]*(double)j/30.+2.-yIntr[point]);
      tempDn->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm_up=(tempUp->ProjectionY("temp",i,i))->Integral();
    double norm_dn=(tempDn->ProjectionY("temp",i,i))->Integral();


    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      tempUp->SetBinContent(i,j,tempUp->GetBinContent(i,j)/norm_up);
      tempDn->SetBinContent(i,j,tempDn->GetBinContent(i,j)/norm_dn);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  histoPair.first  = tempUp;
  histoPair.second = tempDn;

  return histoPair;

}


//=======================================================================
pair<TH2F*,TH2F*> applySystUncForInterference(TH2F* temp){

  // for interference reweighting
  TF1* gauss = new TF1("gauss","gaus",100,1000);

  gauss->SetParameter(0,0.354258);
  gauss->SetParameter(1,114.909);
  gauss->SetParameter(2,17.1512);

  TH2F* tempUp = new TH2F(*temp);
  TH2F* tempDn = new TH2F(*temp);
  
  pair<TH2F*,TH2F*> histoPair(0,0);

  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;

  double slope;

  for(int i=1; i<=temp->GetNbinsX(); i++){

    // choose correct scale factor
    if(i<8){
      slope=.354;
    }else{
      slope=gauss->Eval((double)((i-1)*2+101));
    }

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
      tempUp->SetBinContent(i,j,newTempValue);
      newTempValue = oldTempValue*(1-slope*((double)j/30.-.5));
      tempDn->SetBinContent(i,j,newTempValue);

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm_up=(tempUp->ProjectionY("temp",i,i))->Integral();
    double norm_dn=(tempDn->ProjectionY("temp",i,i))->Integral();

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      tempUp->SetBinContent(i,j,tempUp->GetBinContent(i,j)/norm_up);
      tempDn->SetBinContent(i,j,tempDn->GetBinContent(i,j)/norm_dn);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  histoPair.first  = tempUp;
  histoPair.second = tempDn;

  return histoPair;

}

//=======================================================================
TH2F* reweightForInterference(TH2F* temp){

  // for interference reweighting
  TF1* gauss = new TF1("gauss","gaus",100,1000);

  gauss->SetParameter(0,0.354258);
  gauss->SetParameter(1,114.909);
  gauss->SetParameter(2,17.1512);

  TH2F* newTemp = new TH2F(*temp);
  
  // ---------------------
  // functions for scaling
  // ---------------------
  
  double oldTempValue=0;
  double newTempValue=0;

  double slope;

  for(int i=1; i<=temp->GetNbinsX(); i++){

    // choose correct scale factor
    if(i<8){
      slope=.354;
    }else{
      slope=gauss->Eval((double)((i-1)*2+101));
    }
    
    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      oldTempValue = temp->GetBinContent(i,j);
      newTempValue = oldTempValue*(1+slope*((double)j/30.-.5));
      newTemp->SetBinContent(i,j,newTempValue);


      /*
      if(i<20){
	cout << "---------------------------------" << endl;
	cout << "scaling factor: " << (1+slope*((double)j/30.-.5)) << endl;
	cout << "old value: "  << oldTempValue << endl;
	cout << "new value: "  << newTempValue << endl;
      }
      */

    }// end loop over Y bins

    // -------------- normalize mZZ slice ----------------

    double norm=(newTemp->ProjectionY("temp",i,i))->Integral();

    for(int j=1; j<=temp->GetNbinsY(); j++){
      
      newTemp->SetBinContent(i,j,newTemp->GetBinContent(i,j)/norm);

    }

    // ---------------------------------------------------

  }// end loop over X bins

  return newTemp;

}

//=======================================================================
// builds templates taking fileName from user (wild cards allowed)
//
// channel=4mu,4e,2e2mu
// sample=H*, ZZTo*, ggZZ*
//=======================================================================

void storeLDDistribution(bool signal,char* fileName, char* tag,bool smooth=true){

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    
  vector<TH1F*> vh_LD;

  TFile *file;
  char temp[50];

  TChain* chain = new TChain("SelectedTree");
  chain->Add(fileName);

  if(signal){
    vh_LD = LDDistributionSignal(chain);
    sprintf(temp,"../datafiles/Dsignal_%s.root",tag);
    file= new TFile(temp,"recreate");
  }else{
    vh_LD = LDDistributionBackground(chain);
    sprintf(temp,"../datafiles/Dbackground_%s.root",tag);
    file= new TFile(temp,"recreate");
  }

  cout<<"LD computed. Vector size "<<vh_LD.size()<<endl;

  // using for variable binning
  //TH2F* h_numEvents = new TH2F("numEvents","numEvents",mZZbins-1,binning,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());
  // using for constant binning
  TH2F* h_numEvents = new TH2F("numEvents","numEvents",mZZbins,lowMzz,highMzz,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());

  for (int i=1; i<mZZbins; i++){
    for(int j=1; j<=vh_LD[0]->GetNbinsX(); j++){
      //cout << vh_LD[i-1]->GetBinContent(j) << endl;
      h_numEvents->SetBinContent(i,j,vh_LD[i-1]->GetBinContent(j));
      h_numEvents->SetBinError(i,j,vh_LD[i-1]->GetBinError(j));
    }
  }

  for (int i=1; i<(mZZbins+1); i++){ //the last one is integrated over mzz
    if(vh_LD[i-1]->Integral()>0)
      vh_LD[i-1]->Scale(1./vh_LD[i-1]->Integral());
  }

  // using for variable binning
  //TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",mZZbins-1,binning,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());
  // using for constant binning
  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",mZZbins,lowMzz,highMzz,vh_LD[0]->GetNbinsX(),vh_LD[0]->GetXaxis()->GetXmin(),vh_LD[0]->GetXaxis()->GetXmax());

  for (int i=1; i<mZZbins; i++){
    for(int j=1; j<=vh_LD[0]->GetNbinsX(); j++){
      //cout << vh_LD[i-1]->GetBinContent(j) << endl;
      h_mzzD->SetBinContent(i,j,vh_LD[i-1]->GetBinContent(j));
    }
  }

  TH2F* oldTemp = new TH2F(*h_mzzD);

  if(smooth)
    h_mzzD = smoothTemplate(h_mzzD,h_numEvents);

  pair<TH2F*,TH2F*> histoPair;
  if(!signal)
    histoPair = reweightForCRunc(h_mzzD);
  else{
    if(strcmp(tag,"2e2mu"))
      h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference
    histoPair = applySystUncForInterference(h_mzzD);   // apply systematic unc for lepton interference
  }

  file->cd();
  h_mzzD->Write();
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  file->Close();

}

//=======================================================================
// builds templates similar to storeLDDistribution() but merges 3 final 
// states above ZZ threshold
// channel=4mu,4e,2e2mu
// sample=H*, ZZTo*, ggZZ*
//=======================================================================

void storeLDDistributionV2(char* channel="4mu",int sampleIndex=1, char* dir="/tmp/whitbeck/", bool smooth=true){

  string sample[3]={"H*","ZZTo*","ggZZ*"};
  string tag[3]={"","qqZZ","ggZZ"};

  TChain* singalTree = new TChain("SelectedTree");
  TChain* combTree = new TChain("SelectedTree");
		    
  char singalFile[100];
  char combFile[100];

  sprintf(singalFile,"%s/7TeV_FSR/HZZ%sTree_%s.root",dir,channel,sample[sampleIndex].c_str());
  singalTree->Add(singalFile);
  sprintf(singalFile,"%s/8TeV_FSR/HZZ%sTree_%s.root",dir,channel,sample[sampleIndex].c_str());
  singalTree->Add(singalFile);

  sprintf(combFile,"%s/7TeV_FSR/HZZ*Tree_%s.root",dir,sample[sampleIndex].c_str());
  combTree->Add(combFile);
  sprintf(combFile,"%s/8TeV_FSR/HZZ*Tree_%s.root",dir,sample[sampleIndex].c_str());
  combTree->Add(combFile);
    
  cout << singalFile << endl;
  cout << combFile << endl;

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    
  vector<TH1F*> vh_LD_lowmass;
  vector<TH1F*> vh_LD_highmass;

  TFile *file;
  char temp[50];

  if(sampleIndex==0){
    vh_LD_lowmass  = LDDistributionSignal(singalTree);
    vh_LD_highmass = LDDistributionSignal(combTree);
    sprintf(temp,"../datafiles/Dsignal_%s.root",channel);
    file= new TFile(temp,"recreate");
  }else{
    vh_LD_lowmass  = LDDistributionBackground(singalTree);
    vh_LD_highmass = LDDistributionBackground(combTree);
    sprintf(temp,"../datafiles/Dbackground_%s_%s.root",tag[sampleIndex].c_str(),channel);
    file= new TFile(temp,"recreate");
  }

  cout << "got template slices" << endl;

  //h_numEvents is used for smoothing
  // using for variable binning
  //TH2F* h_numEvents = new TH2F("numEvents","numEvents",mZZbins-1,binning,vh_LD_lowmass[0]->GetNbinsX(),vh_LD_lowmass[0]->GetXaxis()->GetXmin(),vh_LD_lowmass[0]->GetXaxis()->GetXmax());
  // using for constant binning
  TH2F* h_numEvents = new TH2F("numEvents","numEvents",mZZbins,lowMzz,highMzz,vh_LD_lowmass[0]->GetNbinsX(),vh_LD_lowmass[0]->GetXaxis()->GetXmin(),vh_LD_lowmass[0]->GetXaxis()->GetXmax());

  //fill numEvents for low mass
  for (int i=1; i<41; i++){
    for(int j=1; j<=vh_LD_lowmass[0]->GetNbinsX(); j++){
      h_numEvents->SetBinContent(i,j,vh_LD_lowmass[i-1]->GetBinContent(j));
      h_numEvents->SetBinError(i,j,vh_LD_lowmass[i-1]->GetBinError(j));
    }
  }
  //fill numEvents for high mass
  for (int i=41; i<mZZbins; i++){
    for(int j=1; j<=vh_LD_highmass[0]->GetNbinsX(); j++){
      h_numEvents->SetBinContent(i,j,vh_LD_highmass[i-1]->GetBinContent(j));
      h_numEvents->SetBinError(i,j,vh_LD_highmass[i-1]->GetBinError(j));
    }
  }

  cout << "filled h_numEvents" << endl;

  //normalize mzz slices - lowmass
  for (int i=1; i<(mZZbins+1); i++){ //the last one is integrated over mzz
    if(vh_LD_lowmass[i-1]->Integral()>0)
      vh_LD_lowmass[i-1]->Scale(1./vh_LD_lowmass[i-1]->Integral());
  }

  //normalize mzz slices - highmass
  for (int i=1; i<(mZZbins+1); i++){ //the last one is integrated over mzz
    if(vh_LD_highmass[i-1]->Integral()>0)
      vh_LD_highmass[i-1]->Scale(1./vh_LD_highmass[i-1]->Integral());
  }

  cout << "normalized" << endl;

  //using for variable binning
  //TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",mZZbins-1,binning,vh_LD_lowmass[0]->GetNbinsX(),vh_LD_lowmass[0]->GetXaxis()->GetXmin(),vh_LD_lowmass[0]->GetXaxis()->GetXmax());
  // using for constant binning
  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",mZZbins,lowMzz,highMzz,vh_LD_lowmass[0]->GetNbinsX(),vh_LD_lowmass[0]->GetXaxis()->GetXmin(),vh_LD_lowmass[0]->GetXaxis()->GetXmax());

  //fill 2D template for low mass region
  for (int i=1; i<41; i++){
    for(int j=1; j<=vh_LD_lowmass[0]->GetNbinsX(); j++){
      h_mzzD->SetBinContent(i,j,vh_LD_lowmass[i-1]->GetBinContent(j));
    }
  }

  //fill 2D template for high mass region
  for (int i=41; i<mZZbins; i++){
    for(int j=1; j<=vh_LD_highmass[0]->GetNbinsX(); j++){
      h_mzzD->SetBinContent(i,j,vh_LD_highmass[i-1]->GetBinContent(j));
    }
  }

  cout << "filled h_mzzD" << endl;

  TH2F* oldTemp = new TH2F(*h_mzzD);

  cout << "smoothing" << endl;

  if(smooth)
    h_mzzD = smoothTemplate(h_mzzD,h_numEvents);

  pair<TH2F*,TH2F*> histoPair;
  if(sampleIndex>0){
    cout << "adding background syst" << endl;
    histoPair = reweightForCRunc(h_mzzD);
  }else{
    cout << "correcting for interference and adding syst" << endl;
    if(strcmp(channel,"2e2mu"))
      h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference
    histoPair = applySystUncForInterference(h_mzzD);   // apply systematic unc for lepton interference
  }

  file->cd();
  h_mzzD->Write();
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  file->Close();

}

//=======================================================================

TH2F* fillTemplate(char* channel="4mu", int sampleIndex=0,bool isLowMass=true){

  string sample[3]={"H*","ZZTo*","ggZZ*"};
  string sampleName[3]={"signal","qqZZ","ggZZ"};


  TChain* bkgMC = new TChain("SelectedTree");
  char temp[100];
  if(isLowMass){
    sprintf(temp,"CJLSTtrees_June10_2012/7plus8TeV_FSR/HZZ%sTree_%s.root",channel,sample[sampleIndex].c_str());
    bkgMC->Add(temp);
  }else{
    sprintf(temp,"CJLSTtrees_June10_2012/7plus8TeV_FSR/HZZ*Tree_%s.root",sample[sampleIndex].c_str());
  }

  bkgMC->Add(temp);

  float mzz,D,w;

  bkgMC->SetBranchAddress("ZZMass",&mzz);
  bkgMC->SetBranchAddress("ZZLD",&D);
  bkgMC->SetBranchAddress("MC_weight",&w);

  TH2F* bkgHist;
  if(!isLowMass)
    bkgHist = new TH2F("bkgHisto","bkgHisto",310,180,800,30,0,1);
  else
    bkgHist = new TH2F("bkgHisto","bkgHisto",40,100,180,30,0,1);

  bkgHist->Sumw2();

  // fill histogram

  for(int i=0; i<bkgMC->GetEntries(); i++){

    bkgMC->GetEntry(i);

    if(w<.0015)
      bkgHist->Fill(mzz,D,w);

  }

  // normalize slices

  double norm;
  TH1F* tempProj;

  for(int i=1; i<=350; i++){

    tempProj = (TH1F*) bkgHist->ProjectionY("tempProj",i,i);
    norm=tempProj->Integral();

    for(int j=1; j<=30; j++){
      bkgHist->SetBinContent(i,j, bkgHist->GetBinContent(i,j)/norm   );
    }

  }
  
  // average 

  TH2F* notSmooth = new TH2F(*bkgHist);

  if(!isLowMass){
    
    int nXbins=(isLowMass)?40:310;
    int nYbins=30;
    int binMzz;
    int effectiveArea=1;
    
    double average=0,binsUsed=0;
    
    for(int i=1; i<=nXbins; i++){
      for(int j=1; j<=nYbins; j++){
	
	binMzz=(i-1)*2+181;

	if( binMzz<300 ) continue;
	if( binMzz>=300 && binMzz<350 ) effectiveArea=1;
	if( binMzz>=350 && binMzz<500 ) effectiveArea=3;
	if( binMzz>=500 && binMzz<600 ) effectiveArea=5;
	if( binMzz>=600 ) effectiveArea=7;
	
	for(int a=-effectiveArea; a<=effectiveArea; a++){
	  if(a+i<1 || a+i>310 || j>30 || j<1) continue;
	  average+= notSmooth->GetBinContent(a+i,j);
	  binsUsed++;
	}
	
	bkgHist->SetBinContent(i,j,average/binsUsed);

	average=0;
	binsUsed=0;
	
      } // end loop over D
    } // end loop over mZZ
  } // end of horizontal averaging

  // smooth

  bkgHist->Smooth();
  if(!isLowMass)
    bkgHist->Smooth();

  // draw

  TCanvas* can = new TCanvas("can","can",800,400);
  can->Divide(2,1);
  
  can->cd(1);
  bkgHist->Draw("COLZ");
  can->cd(2);
  notSmooth->Draw("COLZ");
  
  if(isLowMass)
    sprintf(temp,"MELAtemplates_smooth_vs_nonSmooth_%s_%s_lowmass.eps",sampleName[sampleIndex].c_str(),channel);
  else
    sprintf(temp,"MELAtemplates_smooth_vs_nonSmooth_%s_%s_highmass.eps",sampleName[sampleIndex].c_str(),channel);

  can->SaveAs(temp);

  if(isLowMass)
    sprintf(temp,"MELAtemplates_smooth_vs_nonSmooth_%s_%s_lowmass.png",sampleName[sampleIndex].c_str(),channel);
  else
    sprintf(temp,"MELAtemplates_smooth_vs_nonSmooth_%s_%s_highmass.png",sampleName[sampleIndex].c_str(),channel);

  can->SaveAs(temp);

  return bkgHist;
  
}

//=======================================================================

TH2F* mergeTemplates(TH2F* lowTemp, TH2F* highTemp){
  
  TH2F* h_mzzD = new TH2F("h_mzzD","h_mzzD",350,100,800,30,0,1);
  
  // copy lowmass into h_mzzD
  for(int i=1; i<=40; i++){
    for(int j=1; j<=30; j++){
      h_mzzD->SetBinContent(i,j, lowTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  // copy high mass into h_mzzD
  for(int i=1; i<=310; i++){
    for(int j=1; j<=30; j++){
      h_mzzD->SetBinContent(i+40,j, highTemp->GetBinContent(i,j)  );
    }// end loop over D
  }// end loop over mZZ

  return h_mzzD;

}

//=======================================================================

void makeTemplate(char* channel="4mu"){

  char temp[150];

  sprintf(temp,"Dsignal_%s.root",channel);
  TFile* fsig = new TFile(temp,"RECREATE");
  sprintf(temp,"Dbackground_qqZZ_%s.root",channel);
  TFile* fqqZZ = new TFile(temp,"RECREATE");
  sprintf(temp,"Dbackground_ggZZ_%s.root",channel);
  TFile* fggZZ = new TFile(temp,"RECREATE");

  TH2F* oldTemp;

  pair<TH2F*,TH2F*> histoPair;

  TH2F* low,*high,*h_mzzD;
  
  low = fillTemplate(channel,0,true);
  high = fillTemplate(channel,0,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "correcting for interference and adding syst" << endl;
  if(strcmp(channel,"2e2mu"))
    h_mzzD = reweightForInterference(h_mzzD);  // correct templates for lepton interference

  cout << "h_mzzD: " << h_mzzD << endl;

  histoPair = applySystUncForInterference(h_mzzD);   // apply systematic unc for lepton interference

  // --------------------------------------------------

  fsig->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fsig->Close();

  // ==========================

  low = fillTemplate(channel,1,true);
  high = fillTemplate(channel,1,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fqqZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fqqZZ->Close();

  // ==========================

  low = fillTemplate(channel,2,true);
  high = fillTemplate(channel,2,false);
  h_mzzD = mergeTemplates(low,high);

  // ---------- apply interference reweighting --------
  
  oldTemp = new TH2F(*h_mzzD);
  oldTemp->SetName("oldTemp");

  cout << "apply systematics for zjets control region" << endl;
  
  histoPair = reweightForCRunc(h_mzzD);

  // --------------------------------------------------

  fggZZ->cd();
  h_mzzD->Write("h_mzzD");
  oldTemp->Write("oldTemp");
  histoPair.first->Write("h_mzzD_up");
  histoPair.second->Write("h_mzzD_dn");
  fggZZ->Close();

}

//=======================================================================

void storeLDDistributionV3(){

  makeTemplate("4mu");
  makeTemplate("4e");
  makeTemplate("2e2mu");

}

//=======================================================================

void makeAllTemplates(){

  storeLDDistribution(true,"../datafiles/HZZ2e2muTree_H*.root","2e2mu");
  storeLDDistribution(true,"../datafiles/HZZ4muTree_H*.root","4mu");
  storeLDDistribution(true,"../datafiles/HZZ4eTree_H*.root","4e");
  
  storeLDDistribution(false,"../datafiles/HZZ2e2muTree_ZZTo*.root","qqZZ_2e2mu");
  storeLDDistribution(false,"../datafiles/HZZ4muTree_ZZTo*.root","qqZZ_4mu");
  storeLDDistribution(false,"../datafiles/HZZ4eTree_ZZTo*.root","qqZZ_4e");
  
  storeLDDistribution(false,"../datafiles/HZZ2e2muTree_ggZZ*.root","ggZZ_2e2mu");
  storeLDDistribution(false,"../datafiles/HZZ4muTree_ggZZ*.root","ggZZ_4mu");
  storeLDDistribution(false,"../datafiles/HZZ4eTree_ggZZ*.root","ggZZ_4e");
  
}

//=======================================================================

void makeTemplatesV2(){

  storeLDDistributionV2("4mu",0);
  storeLDDistributionV2("4mu",1);
  storeLDDistributionV2("4mu",2);

  storeLDDistributionV2("4e",0);
  storeLDDistributionV2("4e",1);
  storeLDDistributionV2("4e",2);

  storeLDDistributionV2("2e2mu",0);
  storeLDDistributionV2("2e2mu",1);
  storeLDDistributionV2("2e2mu",2);

  
}

void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4Lep11, TLorentzVector thep4Lep12, TLorentzVector thep4Z2, TLorentzVector thep4Lep21, TLorentzVector thep4Lep22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2){
	
	
	//std::cout << "In calculate angles..." << std::endl;
	
	double norm;
	
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( thep4Z1 );
	TLorentzVector thep4Z2inXFrame( thep4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
	
	// calculate phi1, phi2, costhetastar
	phi1 = theZ1X_p3.Phi();
	phi2 = theZ2X_p3.Phi();
	
	///////////////////////////////////////////////
	// check for z1/z2 convention, redefine all 4 vectors with convention
	///////////////////////////////////////////////	
	TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;

	/* old convention of choosing Z1 ------------------------------
	p4H = thep4H;
	if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){
		p4Z1 = thep4Z2; p4M11 = thep4Lep21; p4M12 = thep4Lep22;
		p4Z2 = thep4Z1; p4M21 = thep4Lep11; p4M22 = thep4Lep12;		
		costhetastar = theZ2X_p3.CosTheta();
	}
	else{
		p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
		p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
		costhetastar = theZ1X_p3.CosTheta();
	} ---------------------------------------------- */

	p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
	p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
	costhetastar = theZ1X_p3.CosTheta();
	
	//std::cout << "phi1: " << phi1 << ", phi2: " << phi2 << std::endl;
	
	// now helicity angles................................
	// ...................................................
	TVector3 boostZ1 = -(p4Z1.BoostVector());
	TLorentzVector p4Z2Z1(p4Z2);
	p4Z2Z1.Boost(boostZ1);
	//find the decay axis
	/////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
	TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
	norm = 1/(unitx_1.Mag());
	unitx_1*=norm;
	//boost daughters of z2
	TLorentzVector p4M21Z1(p4M21);
	TLorentzVector p4M22Z1(p4M22);
	p4M21Z1.Boost(boostZ1);
	p4M22Z1.Boost(boostZ1);
	//create z and y axes
	/////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
	TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
	TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
	TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
	norm = 1/(unitz_1.Mag());
	unitz_1 *= norm;
	TVector3 unity_1 = unitz_1.Cross(unitx_1);
	
	//caculate theta1
	TLorentzVector p4M11Z1(p4M11);
	p4M11Z1.Boost(boostZ1);
	TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
	TVector3 unitM11 = p3M11.Unit();
	double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
	TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
	costheta1 = M11_Z1frame.CosTheta();
	//std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
	//////-----------------------old way of calculating phi---------------/////////
	phi = M11_Z1frame.Phi();
	
	//set axes for other system
	TVector3 boostZ2 = -(p4Z2.BoostVector());
	TLorentzVector p4Z1Z2(p4Z1);
	p4Z1Z2.Boost(boostZ2);
	TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
	norm = 1/(unitx_2.Mag());
	unitx_2*=norm;
	//boost daughters of z2
	TLorentzVector p4M11Z2(p4M11);
	TLorentzVector p4M12Z2(p4M12);
	p4M11Z2.Boost(boostZ2);
	p4M12Z2.Boost(boostZ2);
	TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
	TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
	TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
	norm = 1/(unitz_2.Mag());
	unitz_2*=norm;
	TVector3 unity_2 = unitz_2.Cross(unitx_2);
	//calcuate theta2
	TLorentzVector p4M21Z2(p4M21);
	p4M21Z2.Boost(boostZ2);
	TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
	TVector3 unitM21 = p3M21.Unit();
	double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
	TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
	costheta2 = M21_Z2frame.CosTheta();
	
	// calculate phi
	//calculating phi_n
	TLorentzVector n_p4Z1inXFrame( p4Z1 );
	TLorentzVector n_p4M11inXFrame( p4M11 );
	n_p4Z1inXFrame.Boost( boostX );
	n_p4M11inXFrame.Boost( boostX );        
	TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
	TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
	TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );
	TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
	TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );
	
	TLorentzVector n_p4M21inXFrame( p4M21 );
	n_p4M21inXFrame.Boost( boostX );
	TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
	//rotate into other plane
	TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
	
	///////-----------------new way of calculating phi-----------------///////
	//double phi_n =  n_p4M21inXFrame_unitprime.Phi();
	/// and then calculate phistar1
	TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
	TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
	// negative sign is for arrow convention in paper
	phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
	
	// and the calculate phistar2
	TLorentzVector n_p4Z2inXFrame( p4Z2 );
	n_p4Z2inXFrame.Boost( boostX );
	TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
	///////TLorentzVector n_p4M21inXFrame( p4M21 );
	//////n_p4M21inXFrame.Boost( boostX );        
	////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
	TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );
	TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
	TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
	TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
	phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());
	
	double phistar12_0 = phistar1 + phistar2;
	if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
	else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
	else phistar12 = phistar12_0;
	
}

