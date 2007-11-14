
#ifndef AnalysisLayers_h
#define AnalysisLayers_h

#include "IOCombined.h"
#include "IOEventTree.h"
#include <vector>

class AnalysisLayers : public IOCombined {

 public:

  // constructor
  //AnalysisLayers(IOCombined* ioc, Bool_t par);
  AnalysisLayers(IOCombined* ioc);

  // destructor
  ~AnalysisLayers();

  // analysis
  void analysis(TFile * rootFile, TCanvas* can, TString psname, TString selstr, TString floatVar, Bool_t PxZoom, 
    Int_t minhit, Float_t* ms, Float_t* mr, Int_t* itp, Int_t nb, Int_t nitermost );
  void globalTitle(TString s);
  void cleanup(void);

 private:

  Int_t nbin,nbin2;
  Double_t maxshift[3],maxshift2[3];
  Double_t maxrot[3],maxrot2[3];
  Float_t maxr,maxz;
  Int_t nlay,laymin,laymax;

  static const Float_t c2m=10000.; // conversion cm->mu

  TObjArray* hstore;

  Bool_t parmode;
};

//-----------------------------------------------------------------------------

AnalysisLayers::AnalysisLayers(IOCombined* ioc) : 
  IOCombined(*ioc) //, 
  //  parmode(par)
{
  hstore= new TObjArray(200,0);
}

//-----------------------------------------------------------------------------

AnalysisLayers::~AnalysisLayers()
{
  cleanup();
}

//-----------------------------------------------------------------------------

void AnalysisLayers::cleanup(void)
{
  for (Int_t i=0;i<1000;i++) {
    if (hstore->At(i)!=0) delete hstore->At(i);
  }
}

//-----------------------------------------------------------------------------

void AnalysisLayers::analysis(TFile * rootFile, TCanvas* can,TString psname, TString selstr, TString floatVar, 
  Bool_t PxZoom, Int_t minhit,  Float_t* ms, Float_t* mr, Int_t* itp, Int_t nb, Int_t nitermost)
{
  
  rootFile->cd("/") ;
  rootFile->mkdir(selstr) ;
  rootFile->cd(selstr) ;

  // iteration numbers used for legend
  //Int_t itpp[4];
  //for (Int_t i=0;i<4;i++) itpp[i]=itp[i];

  //  // adjustments for parallel mode
  //if (parmode) {
  //  for (Int_t i=0;i<4;i++) itp[i]++;
  // }

  // general (number of bins etc)
  nbin=nb; nbin2=nb*2; 
  // nitermost=itp[3]+1;

  if (PxZoom) {
    maxr=15.; maxz=60.; // pixel barrel
    nlay=12; laymin=-6;  laymax=6;      // pixel only
  }
  else {
    maxr=85.; maxz=300.; // full  (originally maxr = 130.)
    nlay=54; laymin=-27; laymax=27;  // full
  }

  // which iterations used for slice plots
  cout <<"Iterations used: "<<itp[0]<<","<<itp[1]<<","<<itp[2]<<","
       <<itp[3]<<endl;

  // max size of shifts
  maxshift[0]=ms[0];  maxshift[1]=ms[0];  maxshift[2]=ms[0];
  maxshift2[0]=ms[1]; maxshift2[1]=ms[1]; maxshift2[2]=ms[1];
  // max size of rotations
  maxrot[0]=mr[0];  maxrot[1]=mr[0];  maxrot[2]=mr[0];
  maxrot2[0]=mr[1];  maxrot2[1]=mr[1];  maxrot2[2]=mr[1];

  can->SetTitle(selstr);

   // selection

   bool sel[ALIMAX];
   for (Int_t ia=0;ia<NAli;ia++) {
     sel[ia]=true;
     if (AliUse[ia]==0) sel[ia]=false;
     if (AliNhit[ia]<minhit) sel[ia]=false;
     if (selstr=="TIB_DS") { // only TIB DL layers
       if (isDoubleSided(ia)==kFALSE) sel[ia]=false;
       if (TMath::Abs(AliType[ia])!=3) sel[ia]=false; 
     }
     else if (selstr=="TOB_DS") { // only TOB DL layers
       if (isDoubleSided(ia)==kFALSE) sel[ia]=false;
       if (TMath::Abs(AliType[ia])!=5) sel[ia]=false; 
     }
     else if (selstr=="TIB_SS") { // only TIB SL layers
       if (isDoubleSided(ia)==kTRUE) sel[ia]=false;
       if (TMath::Abs(AliType[ia])!=3) sel[ia]=false; 
     }
     else if (selstr=="TOB_SS") { // only TOB SL layers
       if (isDoubleSided(ia)==kTRUE) sel[ia]=false;
       if (TMath::Abs(AliType[ia])!=5) sel[ia]=false; 
     }
     else if (selstr=="TIB") { // only TIB
       if (TMath::Abs(AliType[ia])!=3) sel[ia]=false; 
     }
     else if (selstr=="TOB") { // only TOB
       if (TMath::Abs(AliType[ia])!=5) sel[ia]=false; 
     }
     else if (selstr=="TID") { // only TID
       if (TMath::Abs(AliType[ia])!=4) sel[ia]=false; 
     }
     else if (selstr=="TEC") { // only TEC
       if (TMath::Abs(AliType[ia])!=6) sel[ia]=false; 
     }
     else if (selstr=="PXB") { // only PXBarrel
       if (TMath::Abs(AliType[ia])!=1) sel[ia]=false; 
     }
     else if (selstr=="PXEC") { // only PXEndcap
       if (TMath::Abs(AliType[ia])!=2) sel[ia]=false; 
     }
     else if (selstr=="PX") { // only Pixel
       if (TMath::Abs(AliType[ia])>2) sel[ia]=false; 
     }
     else if (selstr=="STRIP") { // only strip
       if (TMath::Abs(AliType[ia])<=2) sel[ia]=false; 
     }
   }


  TString isempty=""; // suppress plotting of global title


  // book histos ------------------------------------------

  // number of hits
  TH1F* hnhit = new TH1F("hnhit"," log_{10}(#Hits); log_{10}Hits",nbin,1.,5.);
  hstore->Add(hnhit);

  TH1F* hstat1=new TH1F("hstat1","",nlay,laymin,laymax);
  TH1F* hstat2=new TH1F("hstat2","",nlay,laymin,laymax);
  TH1F* hstat3=new TH1F("hstat3","Average #Hits per Alignable;Layer",nlay,laymin,laymax);


  //shifts vs iteration
  TH1F *hposdiffx1 = new TH1F("hposdiffx1","; #Deltax [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffx1);
  TH1F *hposdiffx2 = new TH1F("hposdiffx2","; #Deltax [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffx2);
  TH1F *hposdiffx3 = new TH1F("hposdiffx3","; #Deltax [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffx3);
  TH1F *hposdiffx4 = new TH1F("hposdiffx4","; #Deltax [#mum] ; N(Alignables)", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffx4);
  TH1F *hposdiffy1 = new TH1F("hposdiffy1","; #Deltay [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffy1);
  TH1F *hposdiffy2 = new TH1F("hposdiffy2"," ;#Deltay [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffy2);
  TH1F *hposdiffy3 = new TH1F("hposdiffy3"," ;#Deltay [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffy3);
  TH1F *hposdiffy4 = new TH1F("hposdiffy4","; #Deltay [#mum] ; N(Alignables)", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffy4);
  TH1F *hposdiffz1 = new TH1F("hposdiffz1","; #Deltaz [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffz1);
  TH1F *hposdiffz2 = new TH1F("hposdiffz2","; #Deltaz [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffz2);
  TH1F *hposdiffz3 = new TH1F("hposdiffz3","; #Deltaz [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffz3);
  TH1F *hposdiffz4 = new TH1F("hposdiffz4","; #Deltaz [#mum] ; N(Alignables)", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffz4);
  //local shifts vs iteration
  TH1F *hposdiffxl1 = new TH1F("hposdiffxl1","#Deltaxl [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffxl1);
  TH1F *hposdiffxl2 = new TH1F("hposdiffxl2","#Deltaxl [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffxl2);
  TH1F *hposdiffxl3 = new TH1F("hposdiffxl3","#Deltaxl [#mum]", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffxl3);
  TH1F *hposdiffxl4 = new TH1F("hposdiffxl4","#Deltaxl [#mum];#Deltaxl", nbin, -maxshift[0], maxshift[0]);
  hstore->Add(hposdiffxl4);
  TH1F *hposdiffyl1 = new TH1F("hposdiffyl1","#Deltayl [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffyl1);
  TH1F *hposdiffyl2 = new TH1F("hposdiffyl2","#Deltayl [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffyl2);
  TH1F *hposdiffyl3 = new TH1F("hposdiffyl3","#Deltayl [#mum]", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffyl3);
  TH1F *hposdiffyl4 = new TH1F("hposdiffyl4","#Deltayl [#mum];#Deltayl", nbin, -maxshift[1], maxshift[1]);
  hstore->Add(hposdiffyl4);
  TH1F *hposdiffzl1 = new TH1F("hposdiffzl1","#Deltazl [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffzl1);
  TH1F *hposdiffzl2 = new TH1F("hposdiffzl2","#Deltazl [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffzl2);
  TH1F *hposdiffzl3 = new TH1F("hposdiffzl3","#Deltazl [#mum]", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffzl3);
  TH1F *hposdiffzl4 = new TH1F("hposdiffzl4","#Deltazl [#mum];#Deltazl", nbin, -maxshift[2], maxshift[2]);
  hstore->Add(hposdiffzl4);

  // rots vs iteration
  TH1F *hrotdiffa1 = new TH1F("hrotdiffa1","#Deltaa [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffa1);
  TH1F *hrotdiffa2 = new TH1F("hrotdiffa2","#Deltaa [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffa2);
  TH1F *hrotdiffa3 = new TH1F("hrotdiffa3","#Deltaa [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffa3);
  TH1F *hrotdiffa4 = new TH1F("hrotdiffa4","#Deltaa [rad];#Deltaa", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffa4);
  TH1F *hrotdiffb1 = new TH1F("hrotdiffb1","#Deltab [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffb1);
  TH1F *hrotdiffb2 = new TH1F("hrotdiffb2","#Deltab [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffb2);
  TH1F *hrotdiffb3 = new TH1F("hrotdiffb3","#Deltab [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffb3);
  TH1F *hrotdiffb4 = new TH1F("hrotdiffb4","#Deltab [rad];#Deltab", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffb4);
  TH1F *hrotdiffc1 = new TH1F("hrotdiffc1","#Deltac [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffc1);
  TH1F *hrotdiffc2 = new TH1F("hrotdiffc2","#Deltac [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffc2);
  TH1F *hrotdiffc3 = new TH1F("hrotdiffc3","#Deltac [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffc3);
  TH1F *hrotdiffc4 = new TH1F("hrotdiffc4","#Deltac [rad];#Deltac", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffc4);

  // local rots vs iteration
  TH1F *hrotdiffal1 = new TH1F("hrotdiffal1","#Deltaal [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffal1);
  TH1F *hrotdiffal2 = new TH1F("hrotdiffal2","#Deltaal [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffal2);
  TH1F *hrotdiffal3 = new TH1F("hrotdiffal3","#Deltaal [rad]", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffal3);
  TH1F *hrotdiffal4 = new TH1F("hrotdiffal4","#Deltaal [rad];#Deltaal", nbin, -maxrot[0], maxrot[0]);
  hstore->Add(hrotdiffal4);
  TH1F *hrotdiffbl1 = new TH1F("hrotdiffbl1","#Deltabl [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffbl1);
  TH1F *hrotdiffbl2 = new TH1F("hrotdiffbl2","#Deltabl [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffbl2);
  TH1F *hrotdiffbl3 = new TH1F("hrotdiffbl3","#Deltabl [rad]", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffbl3);
  TH1F *hrotdiffbl4 = new TH1F("hrotdiffbl4","#Deltabl [rad];#Deltabl", nbin, -maxrot[1], maxrot[1]);
  hstore->Add(hrotdiffbl4);
  TH1F *hrotdiffcl1 = new TH1F("hrotdiffcl1","#Deltacl [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffcl1);
  TH1F *hrotdiffcl2 = new TH1F("hrotdiffcl2","#Deltacl [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffcl2);
  TH1F *hrotdiffcl3 = new TH1F("hrotdiffcl3","#Deltacl [rad]", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffcl3);
  TH1F *hrotdiffcl4 = new TH1F("hrotdiffcl4","#Deltacl [rad];#Deltacl", nbin, -maxrot[2], maxrot[2]);
  hstore->Add(hrotdiffcl4);

  // shifts vs layer,r,z
  TH2F *hdxvslay = new TH2F("hdxvslay","dx;Layer",nlay,laymin,laymax,nbin2,-maxshift2[0],maxshift2[0]);
  hstore->Add(hdxvslay);
  TH2F *hdyvslay = new TH2F("hdyvslay","dy;Layer",nlay,laymin,laymax,nbin2,-maxshift2[1],maxshift2[1]);
  hstore->Add(hdyvslay);
  TH2F *hdzvslay = new TH2F("hdzvslay","dz;Layer",nlay,laymin,laymax,nbin2,-maxshift2[2],maxshift2[2]);
  hstore->Add(hdzvslay);
  TH2F *hdxvsz = new TH2F("hdxvsz","dx vs z;z [cm];dx",nbin2,-maxz,maxz,nbin2,-maxshift2[0],maxshift2[0]);
  hstore->Add(hdxvsz);
  TH2F *hdyvsz = new TH2F("hdyvsz","dy vs z;z [cm];dy",nbin2,-maxz,maxz,nbin2,-maxshift2[1],maxshift2[1]);
  hstore->Add(hdyvsz);
  TH2F *hdzvsz = new TH2F("hdzvsz","dz vs z;z [cm];dz",nbin2,-maxz,maxz,nbin2,-maxshift2[2],maxshift2[2]);
  hstore->Add(hdzvsz);
  TH2F *hdxvsr = new TH2F("hdxvsr","dx vs R;R [cm];dx",nbin2,0,maxr,nbin2,-maxshift2[0],maxshift2[0]);
  hstore->Add(hdxvsr);
  TH2F *hdyvsr = new TH2F("hdyvsr","dy vs R;R [cm];dy",nbin2,0,maxr,nbin2,-maxshift2[1],maxshift2[1]);
  hstore->Add(hdyvsr);
  TH2F *hdzvsr = new TH2F("hdzvsr","dz vs R;R [cm];dz",nbin2,0,maxr,nbin2,-maxshift2[2],maxshift2[2]);
  // my own addition
  TH2F *hdphivsr = new TH2F("hdphivsr","d#phi vs R;R [cm];dz",nbin2,0,maxr,nbin2,-0.005,0.010);
  TH2F *hdphivsrl1 = new TH2F("hdphivsrl1","d#phi vs R;R [cm];dz",20,22.,30.,nbin2,-0.005,0.010);
  TH2F *hdphivsrl2 = new TH2F("hdphivsrl2","d#phi vs R;R [cm];dz",20,30.,38.,nbin2,-0.005,0.010);
  TH2F *hdphivsrl3 = new TH2F("hdphivsrl3","d#phi vs R;R [cm];dz",20,38.,46.,nbin2,-0.005,0.010);
  TH2F *hdphivsrl4 = new TH2F("hdphivsrl4","d#phi vs R;R [cm];dz",20,46.,56.,nbin2,-0.005,0.010);
  //
  hstore->Add(hdzvsr);
  TH2F *hdxvsph = new TH2F("hdxvsph","dx vs Phi;Phi;dx",nbin2,-4,4,nbin2,-maxshift2[0],maxshift2[0]);
  hstore->Add(hdxvsph);
  TH2F *hdyvsph = new TH2F("hdyvsph","dy vs Phi;Phi;dy",nbin2,-4,4,nbin2,-maxshift2[1],maxshift2[1]);
  hstore->Add(hdyvsph);
  TH2F *hdzvsph = new TH2F("hdzvsph","dz vs Phi;Phi;dz",nbin2,-4,4,nbin2,-maxshift2[2],maxshift2[2]);
  hstore->Add(hdzvsph);
  TH2F *hdxlocvsph = new TH2F("hdxlocvsph","dxloc vs Phi;Phi;dxl",nbin2,0.5,3,nbin2,-maxshift2[2],maxshift2[2]);
  hstore->Add(hdxlocvsph);

  TProfile *hdxvshit = new TProfile("hdxvshit","|dx|;#sqrt{Hits}",10,0.,500.,"");
  hstore->Add(hdxvslay);
  TProfile *hdyvshit = new TProfile("hdyvshit","|dy|;#sqrt{Hits}",10,0.,500.,"");
  hstore->Add(hdyvslay);
  TProfile *hdzvshit = new TProfile("hdzvshit","|dz|;#sqrt{Hits}",10,0.,500.,"");
  hstore->Add(hdzvslay);

  // rots vs layer,r,z
  TH2F *hdavslay = new TH2F("hdavslay","da;Layer;da",nlay,laymin,laymax,nbin2,-maxrot2[0],maxrot2[0]);
  hstore->Add(hdavslay);
  TH2F *hdbvslay = new TH2F("hdbvslay","db;Layer;db",nlay,laymin,laymax,nbin2,-maxrot2[1],maxrot2[1]);
  hstore->Add(hdbvslay);
  TH2F *hdcvslay = new TH2F("hdcvslay","dc;Layer;dc",nlay,laymin,laymax,nbin2,-maxrot2[2],maxrot2[2]);
  hstore->Add(hdcvslay);
  TH2F *hdavsz = new TH2F("hdavsz","da vs z;z;da",nbin2,-maxz,maxz,nbin2,-maxrot2[0],maxrot2[0]);
  hstore->Add(hdavsz);
  TH2F *hdbvsz = new TH2F("hdbvsz","db vs z;z;db",nbin2,-maxz,maxz,nbin2,-maxrot2[1],maxrot2[1]);
  hstore->Add(hdbvsz);
  TH2F *hdcvsz = new TH2F("hdcvsz","dc vs z;z;dc",nbin2,-maxz,maxz,nbin2,-maxrot2[2],maxrot2[2]);
  hstore->Add(hdcvsz);
  TH2F *hdavsr = new TH2F("hdavsr","da vs R;R;dz",nbin2,0,maxr,nbin2,-maxrot2[0],maxrot2[0]);
  hstore->Add(hdavsr);
  TH2F *hdbvsr = new TH2F("hdbvsr","db vs R;R;db",nbin2,0,maxr,nbin2,-maxrot2[1],maxrot2[1]);
  hstore->Add(hdbvsr);
  TH2F *hdcvsr = new TH2F("hdcvsr","dc vs R;R;dc",nbin2,0,maxr,nbin2,-maxrot2[2],maxrot2[2]);
  hstore->Add(hdcvsr);
  TH2F *hdavsph = new TH2F("hdavsph","da vs Phi;Phi;da",nbin2,-4,4,nbin2,-maxrot2[0],maxrot2[0]);
  hstore->Add(hdavsph);
  TH2F *hdbvsph = new TH2F("hdbvsph","db vs Phi;Phi;db",nbin2,-4,4,nbin2,-maxrot2[1],maxrot2[1]);
  hstore->Add(hdbvsph);
  TH2F *hdcvsph = new TH2F("hdcvsph","dc vs Phi;Phi;dc",nbin2,-4,4,nbin2,-maxrot2[2],maxrot2[2]);
  hstore->Add(hdcvsph);

  // Correlation between floated variables
  TH2F *hdxvsdy = new TH2F("hdxvsdy","dy;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxshift[1],maxshift[1]);
  hstore->Add(hdxvsdy);
  TH2F *hdxvsdz = new TH2F("hdxvsdz","dz;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxshift[2],maxshift[2]);
  hstore->Add(hdxvsdz);
  TH2F *hdxvsda = new TH2F("hdxvsda","da;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hdxvsda);
  TH2F *hdxvsdb = new TH2F("hdxvsdb","db;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hdxvsdb);
  TH2F *hdxvsdc = new TH2F("hdxvsdc","dc;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hdxvsdc);
  TH2F *hdyvsdz = new TH2F("hdyvsdz","dz;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxshift[2],maxshift[2]);
  hstore->Add(hdyvsdz);
  TH2F *hdyvsda = new TH2F("hdyvsda","da;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hdyvsda);
  TH2F *hdyvsdb = new TH2F("hdyvsdb","db;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hdyvsdb);
  TH2F *hdyvsdc = new TH2F("hdyvsdc","dc;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hdyvsdc);
  TH2F *hdzvsda = new TH2F("hdzvsda","da;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hdzvsda);
  TH2F *hdzvsdb = new TH2F("hdzvsdb","db;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hdzvsdb);
  TH2F *hdzvsdc = new TH2F("hdzvsdc","dc;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hdzvsdc);
  TH2F *hdavsdb = new TH2F("hdavsdb","db;da",nbin2,-maxrot[0],maxrot[0],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hdavsdb);
  TH2F *hdavsdc = new TH2F("hdavsdc","dc;da",nbin2,-maxrot[0],maxrot[0],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hdavsdc);
  TH2F *hdbvsdc = new TH2F("hdbvsdc","dc;db",nbin2,-maxrot[1],maxrot[1],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hdbvsdc);
  
  // Correlation between steps in floated variables
  TH2F *hstepdxvsdy = new TH2F("hstepdxvsdy","dy;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxshift[1],maxshift[1]);
  hstore->Add(hstepdxvsdy);
  TH2F *hstepdxvsdz = new TH2F("hstepdxvsdz","dz;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxshift[2],maxshift[2]);
  hstore->Add(hstepdxvsdz);
  TH2F *hstepdxvsda = new TH2F("hstepdxvsda","da;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hstepdxvsda);
  TH2F *hstepdxvsdb = new TH2F("hstepdxvsdb","db;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hstepdxvsdb);
  TH2F *hstepdxvsdc = new TH2F("hstepdxvsdc","dc;dx",nbin2,-maxshift[0],maxshift[0],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hstepdxvsdc);
  TH2F *hstepdyvsdz = new TH2F("hstepdyvsdz","dz;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxshift[2],maxshift[2]);
  hstore->Add(hstepdyvsdz);
  TH2F *hstepdyvsda = new TH2F("hstepdyvsda","da;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hstepdyvsda);
  TH2F *hstepdyvsdb = new TH2F("hstepdyvsdb","db;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hstepdyvsdb);
  TH2F *hstepdyvsdc = new TH2F("hstepdyvsdc","dc;dy",nbin2,-maxshift[1],maxshift[1],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hstepdyvsdc);
  TH2F *hstepdzvsda = new TH2F("hstepdzvsda","da;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[0],maxrot[0]);
  hstore->Add(hstepdzvsda);
  TH2F *hstepdzvsdb = new TH2F("hstepdzvsdb","db;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hstepdzvsdb);
  TH2F *hstepdzvsdc = new TH2F("hstepdzvsdc","dc;dz",nbin2,-maxshift[2],maxshift[2],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hstepdzvsdc);
  TH2F *hstepdavsdb = new TH2F("hstepdavsdb","db;da",nbin2,-maxrot[0],maxrot[0],nbin2,-maxrot[1],maxrot[1]);
  hstore->Add(hstepdavsdb);
  TH2F *hstepdavsdc = new TH2F("hstepdavsdc","dc;da",nbin2,-maxrot[0],maxrot[0],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hstepdavsdc);
  TH2F *hstepdbvsdc = new TH2F("hstepdbvsdc","dc;db",nbin2,-maxrot[1],maxrot[1],nbin2,-maxrot[2],maxrot[2]);
  hstore->Add(hstepdbvsdc);

  Float_t maxres=50.;
  Int_t resiter=1;

  TH1F *hresp1 = new TH1F("hresp1","Residual dx",nbin,-maxres,maxres);
  TH1F *hresp2 = new TH1F("hresp2","Residual dy",nbin,-maxres,maxres);
  TH1F *hresp3 = new TH1F("hresp3","Residual dz",nbin,-maxres,maxres);
  TH1F *hresp4 = new TH1F("hresp4","Residual da",nbin,-maxres,maxres);
  TH1F *hresp5 = new TH1F("hresp5","Residual db",nbin,-maxres,maxres);
  TH1F *hresp6 = new TH1F("hresp6","Residual dc",nbin,-maxres,maxres);

  // polyline array for shifts
  TObjArray* thePolylineXArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineYArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineZArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineXlArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineYlArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineZlArr = new TObjArray(ALIMAX,0);

  // polyline array for rots
  TObjArray* thePolylineAArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineBArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineCArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineAlArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineBlArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylineClArr = new TObjArray(ALIMAX,0);

  // polyline array for parameters
  TObjArray* thePolylinePxArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylinePyArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylinePzArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylinePaArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylinePbArr = new TObjArray(ALIMAX,0);
  TObjArray* thePolylinePcArr = new TObjArray(ALIMAX,0);

  Int_t npar=0; /// DEBUG
  Int_t j=0;
  Int_t tempLayer[ALIMAX];

  for (Int_t ia=0;ia<NAli;ia++) {
    if (sel[ia]==1) {

    // Add brutal object selection here if you want to cheat 
    // and show only best results 

    // && AliId[ia] != 369214980 && AliId[ia] != 369264132 && AliId[ia] != 369330436 && AliId[ia] != 369346564) {

     double lognhit = log10(max(0.0001,double(AliNhit[ia])));
     double sqrtnhit = sqrt(max(0.0001,double(AliNhit[ia])));
     hnhit->Fill(lognhit);

       if (NIter>=itp[0]) {
         hposdiffx1->Fill(c2m*(AbsPosA[ia][itp[0]][0]-AbsPosO[ia][0]),1.);
         hposdiffy1->Fill(c2m*(AbsPosA[ia][itp[0]][1]-AbsPosO[ia][1]),1.);
         hposdiffz1->Fill(c2m*(AbsPosA[ia][itp[0]][2]-AbsPosO[ia][2]),1.);
         hrotdiffa1->Fill(    (AbsPosA[ia][itp[0]][3]-AbsPosO[ia][3]),1.);
         hrotdiffb1->Fill(    (AbsPosA[ia][itp[0]][4]-AbsPosO[ia][4]),1.);
         hrotdiffc1->Fill(    (AbsPosA[ia][itp[0]][5]-AbsPosO[ia][5]),1.);
         hposdiffxl1->Fill(c2m*(AbsPosAl[ia][itp[0]][0]-AbsPosOl[ia][0]),1.);
         hposdiffyl1->Fill(c2m*(AbsPosAl[ia][itp[0]][1]-AbsPosOl[ia][1]),1.);
         hposdiffzl1->Fill(c2m*(AbsPosAl[ia][itp[0]][2]-AbsPosOl[ia][2]),1.);
         hrotdiffal1->Fill(    (AbsPosAl[ia][itp[0]][3]-AbsPosOl[ia][3]),1.);
         hrotdiffbl1->Fill(    (AbsPosAl[ia][itp[0]][4]-AbsPosOl[ia][4]),1.);
         hrotdiffcl1->Fill(    (AbsPosAl[ia][itp[0]][5]-AbsPosOl[ia][5]),1.);
       }
       if (NIter>=itp[1]) {
         hposdiffx2->Fill(c2m*(AbsPosA[ia][itp[1]][0]-AbsPosO[ia][0]),1.);
         hposdiffy2->Fill(c2m*(AbsPosA[ia][itp[1]][1]-AbsPosO[ia][1]),1.);
         hposdiffz2->Fill(c2m*(AbsPosA[ia][itp[1]][2]-AbsPosO[ia][2]),1.);
         hrotdiffa2->Fill(    (AbsPosA[ia][itp[1]][3]-AbsPosO[ia][3]),1.);
         hrotdiffb2->Fill(    (AbsPosA[ia][itp[1]][4]-AbsPosO[ia][4]),1.);
         hrotdiffc2->Fill(    (AbsPosA[ia][itp[1]][5]-AbsPosO[ia][5]),1.);
         hposdiffxl2->Fill(c2m*(AbsPosAl[ia][itp[1]][0]-AbsPosOl[ia][0]),1.);
         hposdiffyl2->Fill(c2m*(AbsPosAl[ia][itp[1]][1]-AbsPosOl[ia][1]),1.);
         hposdiffzl2->Fill(c2m*(AbsPosAl[ia][itp[1]][2]-AbsPosOl[ia][2]),1.);
         hrotdiffal2->Fill(    (AbsPosAl[ia][itp[1]][3]-AbsPosOl[ia][3]),1.);
         hrotdiffbl2->Fill(    (AbsPosAl[ia][itp[1]][4]-AbsPosOl[ia][4]),1.);
         hrotdiffcl2->Fill(    (AbsPosAl[ia][itp[1]][5]-AbsPosOl[ia][5]),1.);
       }
       if (NIter>=itp[2]) {
         hposdiffx3->Fill(c2m*(AbsPosA[ia][itp[2]][0]-AbsPosO[ia][0]),1.);
         hposdiffy3->Fill(c2m*(AbsPosA[ia][itp[2]][1]-AbsPosO[ia][1]),1.);
         hposdiffz3->Fill(c2m*(AbsPosA[ia][itp[2]][2]-AbsPosO[ia][2]),1.);
         hrotdiffa3->Fill(    (AbsPosA[ia][itp[2]][3]-AbsPosO[ia][3]),1.);
         hrotdiffb3->Fill(    (AbsPosA[ia][itp[2]][4]-AbsPosO[ia][4]),1.);
         hrotdiffc3->Fill(    (AbsPosA[ia][itp[2]][5]-AbsPosO[ia][5]),1.);
         hposdiffxl3->Fill(c2m*(AbsPosAl[ia][itp[2]][0]-AbsPosOl[ia][0]),1.);
         hposdiffyl3->Fill(c2m*(AbsPosAl[ia][itp[2]][1]-AbsPosOl[ia][1]),1.);
         hposdiffzl3->Fill(c2m*(AbsPosAl[ia][itp[2]][2]-AbsPosOl[ia][2]),1.);
         hrotdiffal3->Fill(    (AbsPosAl[ia][itp[2]][3]-AbsPosOl[ia][3]),1.);
         hrotdiffbl3->Fill(    (AbsPosAl[ia][itp[2]][4]-AbsPosOl[ia][4]),1.);
         hrotdiffcl3->Fill(    (AbsPosAl[ia][itp[2]][5]-AbsPosOl[ia][5]),1.);
       }
       if (NIter>=itp[3]) {
         hposdiffx4->Fill(c2m*(AbsPosA[ia][itp[3]][0]-AbsPosO[ia][0]),1.);
         hposdiffy4->Fill(c2m*(AbsPosA[ia][itp[3]][1]-AbsPosO[ia][1]),1.);
         hposdiffz4->Fill(c2m*(AbsPosA[ia][itp[3]][2]-AbsPosO[ia][2]),1.);
         hrotdiffa4->Fill(    (AbsPosA[ia][itp[3]][3]-AbsPosO[ia][3]),1.);
         hrotdiffb4->Fill(    (AbsPosA[ia][itp[3]][4]-AbsPosO[ia][4]),1.);
         hrotdiffc4->Fill(    (AbsPosA[ia][itp[3]][5]-AbsPosO[ia][5]),1.);
         hposdiffxl4->Fill(c2m*(AbsPosAl[ia][itp[3]][0]-AbsPosOl[ia][0]),1.);
         hposdiffyl4->Fill(c2m*(AbsPosAl[ia][itp[3]][1]-AbsPosOl[ia][1]),1.);
         hposdiffzl4->Fill(c2m*(AbsPosAl[ia][itp[3]][2]-AbsPosOl[ia][2]),1.);
         hrotdiffal4->Fill(    (AbsPosAl[ia][itp[3]][3]-AbsPosOl[ia][3]),1.);
         hrotdiffbl4->Fill(    (AbsPosAl[ia][itp[3]][4]-AbsPosOl[ia][4]),1.);
         hrotdiffcl4->Fill(    (AbsPosAl[ia][itp[3]][5]-AbsPosOl[ia][5]),1.);
       }
      
       // Store correlations
       for (Int_t ii=0 ; ii < NIter ;ii++) {
         hdxvsdy->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosOl[ia][0]),
                       c2m*(AbsPosAl[ia][ii][1]-AbsPosOl[ia][1]));
         hdxvsdz->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosOl[ia][0]),
                       c2m*(AbsPosAl[ia][ii][2]-AbsPosOl[ia][2]));
         hdxvsda->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosOl[ia][0]),
                       AbsPosAl[ia][ii][3]-AbsPosOl[ia][3]);
         hdxvsdb->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosOl[ia][0]),
                       AbsPosAl[ia][ii][4]-AbsPosOl[ia][4]);
         hdxvsdc->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosOl[ia][0]),
                       AbsPosAl[ia][ii][5]-AbsPosOl[ia][5]);
         hdyvsdz->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosOl[ia][1]),
                       c2m*(AbsPosAl[ia][ii][2]-AbsPosOl[ia][2]));
         hdyvsda->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosOl[ia][1]),
                       AbsPosAl[ia][ii][3]-AbsPosOl[ia][3]);
         hdyvsdb->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosOl[ia][1]),
                       AbsPosAl[ia][ii][4]-AbsPosOl[ia][4]);
         hdyvsdc->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosOl[ia][1]),
                       AbsPosAl[ia][ii][5]-AbsPosOl[ia][5]);
	 hdzvsda->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosOl[ia][2]),
                       AbsPosAl[ia][ii][3]-AbsPosOl[ia][3]);
         hdzvsdb->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosOl[ia][2]),
                       AbsPosAl[ia][ii][4]-AbsPosOl[ia][4]);
         hdzvsdc->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosOl[ia][2]),
                       AbsPosAl[ia][ii][5]-AbsPosOl[ia][5]);
         hdavsdb->Fill(AbsPosAl[ia][ii][3]-AbsPosOl[ia][3],
                       AbsPosAl[ia][ii][4]-AbsPosOl[ia][4]);
         hdavsdc->Fill(AbsPosAl[ia][ii][3]-AbsPosOl[ia][3],
                       AbsPosAl[ia][ii][5]-AbsPosOl[ia][5]);
         hdbvsdc->Fill(AbsPosAl[ia][ii][4]-AbsPosOl[ia][4],
                       AbsPosAl[ia][ii][5]-AbsPosOl[ia][5]);
       }
       for (Int_t ii=1 ; ii < NIter ;ii++) {
         hstepdxvsdy->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosAl[ia][ii-1][0]),
                       c2m*(AbsPosAl[ia][ii][1]-AbsPosAl[ia][ii-1][1]));
         hstepdxvsdz->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosAl[ia][ii-1][0]),
                       c2m*(AbsPosAl[ia][ii][2]-AbsPosAl[ia][ii-1][2]));
         hstepdxvsda->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosAl[ia][ii-1][0]),
                       AbsPosAl[ia][ii][3]-AbsPosAl[ia][ii-1][3]);
         hstepdxvsdb->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosAl[ia][ii-1][0]),
                       AbsPosAl[ia][ii][4]-AbsPosAl[ia][ii-1][4]);
         hstepdxvsdc->Fill(c2m*(AbsPosAl[ia][ii][0]-AbsPosAl[ia][ii-1][0]),
                       AbsPosAl[ia][ii][5]-AbsPosAl[ia][ii-1][5]);
         hstepdyvsdz->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosAl[ia][ii-1][1]),
                       c2m*(AbsPosAl[ia][ii][2]-AbsPosAl[ia][ii-1][2]));
         hstepdyvsda->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosAl[ia][ii-1][1]),
                       AbsPosAl[ia][ii][3]-AbsPosAl[ia][ii-1][3]);
         hstepdyvsdb->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosAl[ia][ii-1][1]),
                       AbsPosAl[ia][ii][4]-AbsPosAl[ia][ii-1][4]);
         hstepdyvsdc->Fill(c2m*(AbsPosAl[ia][ii][1]-AbsPosAl[ia][ii-1][1]),
                       AbsPosAl[ia][ii][5]-AbsPosAl[ia][ii-1][5]);
	 hstepdzvsda->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosAl[ia][ii-1][2]),
                       AbsPosAl[ia][ii][3]-AbsPosAl[ia][ii-1][3]);
         hstepdzvsdb->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosAl[ia][ii-1][2]),
                       AbsPosAl[ia][ii][4]-AbsPosAl[ia][ii-1][4]);
         hstepdzvsdc->Fill(c2m*(AbsPosAl[ia][ii][2]-AbsPosAl[ia][ii-1][2]),
                       AbsPosAl[ia][ii][5]-AbsPosAl[ia][ii-1][5]);
         hstepdavsdb->Fill(AbsPosAl[ia][ii][3]-AbsPosAl[ia][ii-1][3],
                       AbsPosAl[ia][ii][4]-AbsPosAl[ia][ii-1][4]);
         hstepdavsdc->Fill(AbsPosAl[ia][ii][3]-AbsPosAl[ia][ii-1][3],
                       AbsPosAl[ia][ii][5]-AbsPosAl[ia][ii-1][5]);
         hstepdbvsdc->Fill(AbsPosAl[ia][ii][4]-AbsPosAl[ia][ii-1][4],
                       AbsPosAl[ia][ii][5]-AbsPosAl[ia][ii-1][5]);
       }

       Float_t dx = c2m*(AbsPosA[ia][NIter][0]-AbsPosO[ia][0]);
       Float_t dxloc = c2m*(AbsPosAl[ia][NIter][0]-AbsPosOl[ia][0]);
       Float_t dy = c2m*(AbsPosA[ia][NIter][1]-AbsPosO[ia][1]);
       Float_t dz = c2m*(AbsPosA[ia][NIter][2]-AbsPosO[ia][2]);
       Float_t da =     (AbsPosA[ia][NIter][3]-AbsPosO[ia][3]);
       Float_t db =     (AbsPosA[ia][NIter][4]-AbsPosO[ia][4]);
       Float_t dc =     (AbsPosA[ia][NIter][5]-AbsPosO[ia][5]);
       Int_t   tl = getlayer(ia);
       Float_t r  = AliPosR[ia];
       Float_t z  = AliPosZ[ia];
       Float_t phi= AliPosPhi[ia];

       cout << AliId[ia] << " " << AliObjId[ia] << " " << c2m*(AbsPosAl[ia][itp[3]][0]- AbsPosOl[ia][0]) << " " << c2m*(AbsPosAl[ia][itp[3]][1]- AbsPosOl[ia][1]) << " " << AliLayer[ia] << " " << r << " " << z << " " << phi << endl;

       // residuals in local frame
       // r = (fitted_shift - true_shift ) / error
       Float_t resl[6];
       for(Int_t i=0;i<6;i++) {
	 resl[i]=0.;
	 if (fabs(AliPar[ia][resiter][i])>0) {
           resl[i] = (AliPar[ia][resiter][i]
                      +(AbsPosAl[ia][resiter-1][i]-AbsPosOl[ia][i])
                     )/AliParErr[ia][resiter][i];
	 }
       }
       if (fabs(resl[0])>0.0000001) hresp1->Fill(resl[0]);
       if (fabs(resl[1])>0.0000001) hresp2->Fill(resl[1]);
       if (fabs(resl[2])>0.0000001) hresp3->Fill(resl[2]);
       if (fabs(resl[3])>0.0000001) hresp4->Fill(resl[3]);
       if (fabs(resl[4])>0.0000001) hresp5->Fill(resl[4]);
       if (fabs(resl[5])>0.0000001) hresp6->Fill(resl[5]);

       hstat1->Fill(tl,AliNhit[ia]);
       hstat2->Fill(tl,1);

       hdxvslay->Fill((float)tl,dx,1.);
       hdyvslay->Fill((float)tl,dy,1.);
       hdzvslay->Fill((float)tl,dz,1.);
       hdavslay->Fill((float)tl,da,1.);
       hdbvslay->Fill((float)tl,db,1.);
       hdcvslay->Fill((float)tl,dc,1.);

       hdxvsz ->Fill(  z,dx,1.);
       hdyvsz ->Fill(  z,dy,1.);
       hdzvsz ->Fill(  z,dz,1.);
       //
       float mydphi = sqrt(pow(dx,2)+pow(dy,2))*(dx/fabs(dx))/(c2m*r);
       hdphivsr ->Fill(  r,mydphi,1.);
       hdphivsrl1 ->Fill(  r,mydphi,1.);
       hdphivsrl2 ->Fill(  r,mydphi,1.);
       hdphivsrl3 ->Fill(  r,mydphi,1.);
       hdphivsrl4 ->Fill(  r,mydphi,1.);
       //
       hdxvsr ->Fill(  r,dx,1.);
       hdyvsr ->Fill(  r,dy,1.);
       hdzvsr ->Fill(  r,dz,1.);
       hdxvsph->Fill(phi,dx,1.);
       hdyvsph->Fill(phi,dy,1.);
       hdzvsph->Fill(phi,dz,1.);
       hdxlocvsph->Fill(phi,dxloc,1.);

       hdavsz ->Fill(  z,da,1.);
       hdbvsz ->Fill(  z,db,1.);
       hdcvsz ->Fill(  z,dc,1.);
       hdavsr ->Fill(  r,da,1.);
       hdbvsr ->Fill(  r,db,1.);
       hdcvsr ->Fill(  r,dc,1.);
       hdavsph->Fill(phi,da,1.);
       hdbvsph->Fill(phi,db,1.);
       hdcvsph->Fill(phi,dc,1.);

//       cout << "sqrtnhit" << " " << sqrtnhit << endl;
//       cout << "abs(dx)" << " " << abs(dx) << endl;

       hdxvshit->Fill(sqrtnhit,fabs(dx),1.);
       hdyvshit->Fill(sqrtnhit,fabs(dy),1.);
       hdzvshit->Fill(sqrtnhit,fabs(dz),1.);


       Float_t xarr[ITERMAX];
       Float_t yxarr[ITERMAX],yyarr[ITERMAX],yzarr[ITERMAX];
       Float_t yxlarr[ITERMAX],yylarr[ITERMAX],yzlarr[ITERMAX];
       Float_t yaarr[ITERMAX],ybarr[ITERMAX],ycarr[ITERMAX];
       Float_t yalarr[ITERMAX],yblarr[ITERMAX],yclarr[ITERMAX];
       Float_t ypxarr[ITERMAX],ypyarr[ITERMAX],ypzarr[ITERMAX];
       Float_t ypaarr[ITERMAX],ypbarr[ITERMAX],ypcarr[ITERMAX];

       for (Int_t i=0;i<=NIter;i++) {
         xarr[i]=(i);
         yxarr[i] =max(-maxshift[0],min(maxshift[0],
                   c2m*(AbsPosA[ia][i][0]-AbsPosO[ia][0])));
         yyarr[i] =max(-maxshift[1],min(maxshift[1],
                   c2m*(AbsPosA[ia][i][1]-AbsPosO[ia][1])));
         yzarr[i] =max(-maxshift[2],min(maxshift[2],
                   c2m*(AbsPosA[ia][i][2]-AbsPosO[ia][2])));
         yaarr[i] =max(-maxrot[0],min(maxrot[0],
                       (AbsPosA[ia][i][3]-AbsPosO[ia][3])));
         ybarr[i] =max(-maxrot[1],min(maxrot[1],
                       (AbsPosA[ia][i][4]-AbsPosO[ia][4])));
         ycarr[i] =max(-maxrot[2],min(maxrot[2],
                       (AbsPosA[ia][i][5]-AbsPosO[ia][5])));
         yxlarr[i]=max(-maxshift[0],min(maxshift[0],
                   c2m*(AbsPosAl[ia][i][0]-AbsPosOl[ia][0])));
         yylarr[i]=max(-maxshift[1],min(maxshift[1],
                   c2m*(AbsPosAl[ia][i][1]-AbsPosOl[ia][1])));
         yzlarr[i]=max(-maxshift[2],min(maxshift[2],
                   c2m*(AbsPosAl[ia][i][2]-AbsPosOl[ia][2])));
         yalarr[i]=max(-maxrot[0],min(maxrot[0],
                       (AbsPosAl[ia][i][3]-AbsPosOl[ia][3])));
         yblarr[i]=max(-maxrot[1],min(maxrot[1],
                       (AbsPosAl[ia][i][4]-AbsPosOl[ia][4])));
         yclarr[i]=max(-maxrot[2],min(maxrot[2],
                       (AbsPosAl[ia][i][5]-AbsPosOl[ia][5])));
	 ypxarr[i]=max(-maxshift[0],min(maxshift[0],c2m*AliPar[ia][i][0]));
	 ypyarr[i]=max(-maxshift[1],min(maxshift[1],c2m*AliPar[ia][i][1]));
	 ypzarr[i]=max(-maxshift[2],min(maxshift[2],c2m*AliPar[ia][i][2]));
	 ypaarr[i]=max(-maxrot[0],min(maxrot[0],AliPar[ia][i][3]));
	 ypbarr[i]=max(-maxrot[1],min(maxrot[1],AliPar[ia][i][4]));
	 ypcarr[i]=max(-maxrot[2],min(maxrot[2],AliPar[ia][i][5]));
	
       }
       TPolyLine* plx = new TPolyLine(NIter+1,xarr,yxarr,"S");
       thePolylineXArr->Add(plx);
       TPolyLine* ply = new TPolyLine(NIter+1,xarr,yyarr,"S");
       thePolylineYArr->Add(ply);
       TPolyLine* plz = new TPolyLine(NIter+1,xarr,yzarr,"S");
       thePolylineZArr->Add(plz);
       TPolyLine* pla = new TPolyLine(NIter+1,xarr,yaarr,"S");
       thePolylineAArr->Add(pla);
       TPolyLine* plb = new TPolyLine(NIter+1,xarr,ybarr,"S");
       thePolylineBArr->Add(plb);
       TPolyLine* plc = new TPolyLine(NIter+1,xarr,ycarr,"S");
       thePolylineCArr->Add(plc);
       TPolyLine* plxl = new TPolyLine(NIter+1,xarr,yxlarr,"S");
       thePolylineXlArr->Add(plxl);
       TPolyLine* plyl = new TPolyLine(NIter+1,xarr,yylarr,"S");
       thePolylineYlArr->Add(plyl);
       TPolyLine* plzl = new TPolyLine(NIter+1,xarr,yzlarr,"S");
       thePolylineZlArr->Add(plzl);
       TPolyLine* plal = new TPolyLine(NIter+1,xarr,yalarr,"S");
       thePolylineAlArr->Add(plal);
       TPolyLine* plbl = new TPolyLine(NIter+1,xarr,yblarr,"S");
       thePolylineBlArr->Add(plbl);
       TPolyLine* plcl = new TPolyLine(NIter+1,xarr,yclarr,"S");
       thePolylineClArr->Add(plcl);

       TPolyLine* plpx = new TPolyLine(NIter+1,xarr,ypxarr,"S");
       thePolylinePxArr->Add(plpx);
       TPolyLine* plpy = new TPolyLine(NIter+1,xarr,ypyarr,"S");
       thePolylinePyArr->Add(plpy);
       TPolyLine* plpz = new TPolyLine(NIter+1,xarr,ypzarr,"S");
       thePolylinePzArr->Add(plpz);
       TPolyLine* plpa = new TPolyLine(NIter+1,xarr,ypaarr,"S");
       thePolylinePaArr->Add(plpa);
       TPolyLine* plpb = new TPolyLine(NIter+1,xarr,ypbarr,"S");
       thePolylinePbArr->Add(plpb);
       TPolyLine* plpc = new TPolyLine(NIter+1,xarr,ypcarr,"S");
       thePolylinePcArr->Add(plpc);
       tempLayer[j] = AliLayer[ia];
       j++;

     }
  }

  hstat3->Divide(hstat1,hstat2);

  cout <<"found alignables with non-zero params: " << npar << endl;

   // -------------------------------------------------------------------------
   // -------------------------------------------------------------------------
   // -------------------------------------------------------------------------

   can->Clear();
   globalTitle(isempty);

   can->Divide(2,2);
   can->cd(1); hnhit->Draw();
   can->cd(2); 
//   gPad->SetLogy(); 
   PlotTLHisto(hstat3); 

   can->Update();
   can->Print(psname,"ps"); 

   // -------------------------------------------------------------------------
   // shifts ------------------------------------------------------------------

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   // legend for slice plots
   TLegend* leg = new TLegend(0.2, 0.5, 0.45, 0.9);
   leg->SetTextSize(0.075);
   leg->SetBorderSize(1);
   leg->SetFillColor(10);
   leg->SetHeader("After I Iterations:");   
   char siter[5];
   sprintf(siter,"I=%d",itp[0]);
   leg->AddEntry(hposdiffx1,siter,"L");
   sprintf(siter,"I=%d",itp[1]);
   leg->AddEntry(hposdiffx2,siter,"L");
   sprintf(siter,"I=%d",itp[2]);
   leg->AddEntry(hposdiffx3,siter,"L");
   sprintf(siter,"I=%d",itp[3]);
   leg->AddEntry(hposdiffx4,siter,"L");

   Int_t nbin2iter=101;

   // now loop again over alignables to fill convergence plots
   TH2F *hposxvsiter =  new TH2F("hposxvsiter","; Iteration ; #Deltax [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[0], maxshift[0]);
   hstore->Add(hposxvsiter);
   TH2F *hposyvsiter =  new TH2F("hposyvsiter","; Iteration ;#Deltay [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[1], maxshift[1]);
   hstore->Add(hposyvsiter);
   TH2F *hposzvsiter =  new TH2F("hposzvsiter","; Iteration ;#Deltaz [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[2], maxshift[2]);
   hstore->Add(hposzvsiter);

   hposxvsiter->SetStats(kFALSE);
   hposyvsiter->SetStats(kFALSE);
   hposzvsiter->SetStats(kFALSE);
   hposdiffx4->SetStats(kTRUE);
   hposdiffy4->SetStats(kTRUE);
   hposdiffz4->SetStats(kTRUE);
   ofstream myfile1;
   char nomefile[30];
   sprintf(nomefile,"meanerr%s.txt",selstr.Data());
   myfile1.open (nomefile);
   
   can->cd(2); 
   hposdiffx4->Draw();
   myfile1 << "x_glob: After alignment  MEAN = " << hposdiffx4->GetMean() << " RMS = " << hposdiffx4->GetRMS() << endl;
   hposdiffx3->SetLineStyle(2); hposdiffx3->SetLineColor(kRed);
   hposdiffx3->Draw("SAME");
   hposdiffx2->SetLineStyle(3); hposdiffx2->SetLineColor(kGreen); 
   hposdiffx2->Draw("SAME"); 
   myfile1 << "        Before alignment MEAN = " << hposdiffx1->GetMean() << " RMS = " << hposdiffx1->GetRMS() << endl;
   hposdiffx1->SetLineStyle(4); hposdiffx1->SetLineColor(kBlue); 
   hposdiffx1->Draw("SAME");
   leg->Draw();

   can->cd(4); 
   hposdiffy4->Draw();
   myfile1 << "y_glob: After alignment  MEAN = " << hposdiffy4->GetMean() << " RMS = " << hposdiffy4->GetRMS() << endl;
   hposdiffy3->SetLineStyle(2); hposdiffy3->SetLineColor(kRed);
   hposdiffy3->Draw("SAME");
   hposdiffy2->SetLineStyle(3); hposdiffy2->SetLineColor(kGreen);
   hposdiffy2->Draw("SAME");
   hposdiffy1->SetLineStyle(4); hposdiffy1->SetLineColor(kBlue);
   hposdiffy1->Draw("SAME");
   myfile1 << "        Before alignment MEAN = " << hposdiffy1->GetMean() << " RMS = " << hposdiffy1->GetRMS() << endl;

   can->cd(6); 
   hposdiffz4->Draw();
   myfile1 << "z_glob: After alignment  MEAN = " << hposdiffz4->GetMean() << " RMS = " << hposdiffz4->GetRMS() << endl;
   hposdiffz3->SetLineStyle(2); hposdiffz3->SetLineColor(kRed);
   hposdiffz3->Draw("SAME");
   hposdiffz2->SetLineStyle(3); hposdiffz2->SetLineColor(kGreen);
   hposdiffz2->Draw("SAME");
   myfile1 << "        Before alignment MEAN = " << hposdiffz1->GetMean() << " RMS = " << hposdiffz1->GetRMS() << endl;
   hposdiffz1->SetLineStyle(4); hposdiffz1->SetLineColor(kBlue);
   hposdiffz1->Draw("SAME");
   // myfile1.close();

   can->cd(1); 
   hposxvsiter->Fill(0.,0.,1.);
   hposxvsiter->Draw("");
   can->cd(3); 
   hposyvsiter->Fill(0.,0.,1.);
   hposyvsiter->Draw("");
   can->cd(5); 
   hposzvsiter->Fill(0.,0.,1.);
   hposzvsiter->Draw("");

    Double_t plotfrac=1.0;
//   Double_t plotfrac=0.2;
   Int_t nplotted=0;
   Int_t lsty=3;
       
    for (Int_t i=0;i<ALIMAX;i++) { 
      TPolyLine* plx = (TPolyLine*) (thePolylineXArr->At(i));
      TPolyLine* ply = (TPolyLine*) (thePolylineYArr->At(i));
      TPolyLine* plz = (TPolyLine*) (thePolylineZArr->At(i));
      lsty = tempLayer[i];
      if (plx!=0 && ply!=0 && plz!=0 && lsty!=0) {
        Double_t rnd = gRandom->Rndm();
        if (rnd<plotfrac) {
	  nplotted++;
          can->cd(1);
          plx->SetLineWidth(1); plx->SetLineStyle(lsty); 
          plx->SetLineColor(kBlack);
          plx->Draw();

          can->cd(3);
          ply->SetLineWidth(1); ply->SetLineStyle(lsty); 
          ply->SetLineColor(kBlack);
          ply->Draw();

          can->cd(5);
          plz->SetLineWidth(1); plz->SetLineStyle(lsty); 
          plz->SetLineColor(kBlack);
          plz->Draw();
        }
      }
    }

    cout <<"plotted: " << nplotted << endl;

   can->Update();
   can->Print(psname,"ps"); 

   // local shifts ------------------------------------------------------------------

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   // now loop again over alignables to fill convergence plots
   TH2F *hposxlvsiter =  new TH2F("hposxlvsiter","Iteration vs #Deltaxl [#mum];Iteration;#Deltaxl [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[0], maxshift[0]);
   hstore->Add(hposxlvsiter);
   TH2F *hposylvsiter =  new TH2F("hposylvsiter","Iteration vs #Deltayl [#mum];Iteration;#Deltayl [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[1], maxshift[1]);
   hstore->Add(hposylvsiter);
   TH2F *hposzlvsiter =  new TH2F("hposzlvsiter","Iteration vs #Deltazl [#mum];Iteration;#Deltazl [#mum]",
     nitermost,0,nitermost,nbin2iter, -maxshift[2], maxshift[2]);
   hstore->Add(hposzlvsiter);

   hposxlvsiter->SetStats(kFALSE);
   hposylvsiter->SetStats(kFALSE);
   hposzlvsiter->SetStats(kFALSE);
   hposdiffxl4->SetStats(kTRUE);
   hposdiffyl4->SetStats(kTRUE);
   hposdiffzl4->SetStats(kTRUE);
   ofstream myfile2;
   sprintf(nomefile,"meanerr_loc%s.txt",selstr.Data());
   myfile2.open (nomefile);

   can->cd(2); 
   hposdiffxl4->Draw();
   myfile2 << "x_loc: After alignment  MEAN = " << hposdiffxl4->GetMean() << " RMS = " << hposdiffxl4->GetRMS() << endl;
   hposdiffxl3->SetLineStyle(2); hposdiffxl3->SetLineColor(kRed);
   hposdiffxl3->Draw("SAME");
   hposdiffxl2->SetLineStyle(3); hposdiffxl2->SetLineColor(kGreen); 
   hposdiffxl2->Draw("SAME");
   hposdiffxl1->SetLineStyle(4); hposdiffxl1->SetLineColor(kBlue); 
   hposdiffxl1->Draw("SAME");
   myfile2 << "x_loc: Before alignment MEAN = " << hposdiffxl1->GetMean() << " RMS = " << hposdiffxl1->GetRMS() << endl;
   leg->Draw();

   can->cd(4); 
   hposdiffyl4->Draw();
   myfile2 << "y_loc: After alignment  MEAN = " << hposdiffyl4->GetMean() << " RMS = " << hposdiffyl4->GetRMS() << endl;
   hposdiffyl3->SetLineStyle(2); hposdiffyl3->SetLineColor(kRed);
   hposdiffyl3->Draw("SAME");
   hposdiffyl2->SetLineStyle(3); hposdiffyl2->SetLineColor(kGreen);
   hposdiffyl2->Draw("SAME");
   hposdiffyl1->SetLineStyle(4); hposdiffyl1->SetLineColor(kBlue);
   hposdiffyl1->Draw("SAME");
   myfile2 << "y_loc: Before alignment MEAN = " << hposdiffyl1->GetMean() << " RMS = " << hposdiffyl1->GetRMS() << endl;

   can->cd(6); 
   hposdiffzl4->Draw();
   myfile2 << "z_loc: After alignment  MEAN = " << hposdiffzl4->GetMean() << " RMS = " << hposdiffzl4->GetRMS() << endl;
   hposdiffzl3->SetLineStyle(2); hposdiffzl3->SetLineColor(kRed);
   hposdiffzl3->Draw("SAME");
   hposdiffzl2->SetLineStyle(3); hposdiffzl2->SetLineColor(kGreen);
   hposdiffzl2->Draw("SAME");
   hposdiffzl1->SetLineStyle(4); hposdiffzl1->SetLineColor(kBlue);
   hposdiffzl1->Draw("SAME");
   myfile2 << "z_loc: Before alignment MEAN = " << hposdiffzl1->GetMean() << " RMS = " << hposdiffzl1->GetRMS() << endl;

   can->cd(1); 
   hposxlvsiter->Fill(0.,0.,1.);
   hposxlvsiter->Draw("");
   can->cd(3); 
   hposylvsiter->Fill(0.,0.,1.);
   hposylvsiter->Draw("");
   can->cd(5); 
   hposzlvsiter->Fill(0.,0.,1.);
   hposzlvsiter->Draw("");
       
    for (Int_t i=0;i<ALIMAX;i++) { 
      can->cd(1);
      lsty = tempLayer[i]; 
      TPolyLine* plx = (TPolyLine*) (thePolylineXlArr->At(i));
      if (plx!=0) { plx->SetLineWidth(1); plx->SetLineStyle(lsty); plx->Draw();}
      can->cd(3);
      TPolyLine* ply = (TPolyLine*) (thePolylineYlArr->At(i));
      if (ply!=0) { ply->SetLineWidth(1); ply->SetLineStyle(lsty); ply->Draw();}
      can->cd(5);
      TPolyLine* plz = (TPolyLine*) (thePolylineZlArr->At(i));
      if (plz!=0) { plz->SetLineWidth(1); plz->SetLineStyle(lsty); plz->Draw();}
    }

   can->Update();
   can->Print(psname,"ps"); 

   // rotations -------------------------------------------------------------------------
   
   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   // now loop again over alignables to fill convergence plots
   //   Int_t nbin2iter=101;
//   TH2F *hrotavsiter =  new TH2F("hrotavsiter","Iteration vs #Deltaa [rad];Iteration;da",
//     nitermost,0,nitermost,nbin2iter, -maxrot[0], maxrot[0]);
//   hstore->Add(hrotavsiter);
//   TH2F *hrotbvsiter =  new TH2F("hrotbvsiter","Iteration vs #Deltab [rad];Iteration;db",
//     nitermost,0,nitermost,nbin2iter, -maxrot[1], maxrot[1]);
//   hstore->Add(hrotbvsiter);
//   TH2F *hrotcvsiter =  new TH2F("hrotcvsiter","Iteration vs #Deltac [rad];Iteration;dc",
//     nitermost,0,nitermost,nbin2iter, -maxrot[2], maxrot[2]);
   TH2F *hrotavsiter =  new TH2F("hrotavsiter",";Iteration;da",
     nitermost,0,nitermost,nbin2iter, -maxrot[0], maxrot[0]);
   hstore->Add(hrotavsiter);
   TH2F *hrotbvsiter =  new TH2F("hrotbvsiter",";Iteration;db",
     nitermost,0,nitermost,nbin2iter, -maxrot[1], maxrot[1]);
   hstore->Add(hrotbvsiter);
   TH2F *hrotcvsiter =  new TH2F("hrotcvsiter",";Iteration;dc",
     nitermost,0,nitermost,nbin2iter, -maxrot[2], maxrot[2]);
   hstore->Add(hrotcvsiter);

   hrotavsiter->SetStats(kFALSE);
   hrotbvsiter->SetStats(kFALSE);
   hrotcvsiter->SetStats(kFALSE);
   hrotdiffa4->SetStats(kTRUE);
   hrotdiffb4->SetStats(kTRUE);
   hrotdiffc4->SetStats(kTRUE);

   can->cd(2); 
   hrotdiffa4->Draw();
   myfile1 << "alpha_glob: After alignment  MEAN = " << hrotdiffa4->GetMean() << " RMS = " << hrotdiffa4->GetRMS() << endl;
   hrotdiffa3->SetLineStyle(2); hrotdiffa3->SetLineColor(kRed);
   hrotdiffa3->Draw("SAME");
   hrotdiffa2->SetLineStyle(3); hrotdiffa2->SetLineColor(kGreen); 
   hrotdiffa2->Draw("SAME");
   hrotdiffa1->SetLineStyle(4); hrotdiffa1->SetLineColor(kBlue); 
   hrotdiffa1->Draw("SAME");
   myfile1 << "            Before alignment MEAN = " << hrotdiffa1->GetMean() << " RMS = " << hrotdiffa1->GetRMS() << endl;
   leg->Draw();

   can->cd(4); 
   myfile1 << "beta_glob:  After alignment  MEAN = " << hrotdiffb4->GetMean() << " RMS = " << hrotdiffb4->GetRMS() << endl; 
   hrotdiffb4->Draw();
   hrotdiffb3->SetLineStyle(2); hrotdiffb3->SetLineColor(kRed);
   hrotdiffb3->Draw("SAME");
   hrotdiffb2->SetLineStyle(3); hrotdiffb2->SetLineColor(kGreen);
   hrotdiffb2->Draw("SAME");
   hrotdiffb1->SetLineStyle(4); hrotdiffb1->SetLineColor(kBlue);
   hrotdiffb1->Draw("SAME");
   myfile1 << "            Before alignment MEAN = " << hrotdiffb1->GetMean() << " RMS = " << hrotdiffb1->GetRMS() << endl;

   can->cd(6); 
   hrotdiffc4->Draw();
   myfile1 << "gamma_glob: After alignment  MEAN = " << hrotdiffc4->GetMean() << " RMS = " << hrotdiffc4->GetRMS() << endl;
   hrotdiffc3->SetLineStyle(2); hrotdiffc3->SetLineColor(kRed);
   hrotdiffc3->Draw("SAME");
   hrotdiffc2->SetLineStyle(3); hrotdiffc2->SetLineColor(kGreen);
   hrotdiffc2->Draw("SAME");
   hrotdiffc1->SetLineStyle(4); hrotdiffc1->SetLineColor(kBlue);
   myfile1 << "            Before alignment MEAN = " << hrotdiffc1->GetMean() << " RMS = " << hrotdiffc1->GetRMS() << endl << endl;
   hrotdiffc1->Draw("SAME");
   myfile1.close();

   can->cd(1); 
   hrotavsiter->Fill(0.,0.,1.);
   hrotavsiter->Draw("");
   can->cd(3); 
   hrotbvsiter->Fill(0.,0.,1.);
   hrotbvsiter->Draw("");
   can->cd(5); 
   hrotcvsiter->Fill(0.,0.,1.);
   hrotcvsiter->Draw("");

   for (Int_t i=0;i<ALIMAX;i++) { 
     can->cd(1);
     lsty = tempLayer[i];
     TPolyLine* plx = (TPolyLine*) (thePolylineAArr->At(i));
     if (plx!=0) { plx->SetLineWidth(1); plx->SetLineStyle(lsty); plx->Draw();}
     can->cd(3);
     TPolyLine* ply = (TPolyLine*) (thePolylineBArr->At(i));
     if (ply!=0) { ply->SetLineWidth(1); ply->SetLineStyle(lsty); ply->Draw();}
     can->cd(5);
     TPolyLine* plz = (TPolyLine*) (thePolylineCArr->At(i));
     if (plz!=0) { plz->SetLineWidth(1); plz->SetLineStyle(lsty); plz->Draw();}
   }

   can->Update();
   can->Print(psname,"ps"); 

   // LOCAL rotations -------------------------------------------------------------------------
   
   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   // now loop again over alignables to fill convergence plots
   TH2F *hrotalvsiter =  new TH2F("hrotalvsiter","Iteration vs #Deltaal [rad];Iteration;dal",
     nitermost,0,nitermost,nbin2iter, -maxrot[0], maxrot[0]);
   hstore->Add(hrotavsiter);
   TH2F *hrotblvsiter =  new TH2F("hrotblvsiter","Iteration vs #Deltabl [rad];Iteration;dbl",
     nitermost,0,nitermost,nbin2iter, -maxrot[1], maxrot[1]);
   hstore->Add(hrotbvsiter);
   TH2F *hrotclvsiter =  new TH2F("hrotclvsiter","Iteration vs #Deltacl [rad];Iteration;dcl",
     nitermost,0,nitermost,nbin2iter, -maxrot[2], maxrot[2]);
   hstore->Add(hrotcvsiter);

   hrotalvsiter->SetStats(kFALSE);
   hrotblvsiter->SetStats(kFALSE);
   hrotclvsiter->SetStats(kFALSE);
   hrotdiffal4->SetStats(kTRUE);
   hrotdiffbl4->SetStats(kTRUE);
   hrotdiffcl4->SetStats(kTRUE);

   can->cd(2); 
   hrotdiffal4->Draw();
   hrotdiffal3->SetLineStyle(2); hrotdiffal3->SetLineColor(kRed);
   hrotdiffal3->Draw("SAME");
   hrotdiffal2->SetLineStyle(3); hrotdiffal2->SetLineColor(kGreen); 
   hrotdiffal2->Draw("SAME");
   hrotdiffal1->SetLineStyle(4); hrotdiffal1->SetLineColor(kBlue); 
   hrotdiffal1->Draw("SAME");
   leg->Draw();

   can->cd(4); 
   hrotdiffbl4->Draw();
   hrotdiffbl3->SetLineStyle(2); hrotdiffbl3->SetLineColor(kRed);
   hrotdiffbl3->Draw("SAME");
   hrotdiffbl2->SetLineStyle(3); hrotdiffbl2->SetLineColor(kGreen);
   hrotdiffbl2->Draw("SAME");
   hrotdiffbl1->SetLineStyle(4); hrotdiffbl1->SetLineColor(kBlue);
   hrotdiffbl1->Draw("SAME");

   can->cd(6); 
   hrotdiffcl4->Draw();
   hrotdiffcl3->SetLineStyle(2); hrotdiffcl3->SetLineColor(kRed);
   hrotdiffcl3->Draw("SAME");
   hrotdiffcl2->SetLineStyle(3); hrotdiffcl2->SetLineColor(kGreen);
   hrotdiffcl2->Draw("SAME");
   hrotdiffcl1->SetLineStyle(4); hrotdiffcl1->SetLineColor(kBlue);
   hrotdiffcl1->Draw("SAME");

   can->cd(1); 
   hrotalvsiter->Fill(0.,0.,1.);
   hrotalvsiter->Draw("");
   can->cd(3); 
   hrotblvsiter->Fill(0.,0.,1.);
   hrotblvsiter->Draw("");
   can->cd(5); 
   hrotclvsiter->Fill(0.,0.,1.);
   hrotclvsiter->Draw("");

   for (Int_t i=0;i<ALIMAX;i++) { 
     can->cd(1);
     lsty = tempLayer[i];
     TPolyLine* plx = (TPolyLine*) (thePolylineAlArr->At(i));
     if (plx!=0) { plx->SetLineWidth(1); plx->SetLineStyle(lsty); plx->Draw();}
     can->cd(3);
     TPolyLine* ply = (TPolyLine*) (thePolylineBlArr->At(i));
     if (ply!=0) { ply->SetLineWidth(1); ply->SetLineStyle(lsty); ply->Draw();}
     can->cd(5);
     TPolyLine* plz = (TPolyLine*) (thePolylineClArr->At(i));
     if (plz!=0) { plz->SetLineWidth(1); plz->SetLineStyle(lsty); plz->Draw();}
   }

   can->Update();
   can->Print(psname,"ps"); 

   // parameters --------------------------------------------------------------

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   //maxshift2[0]=1;
   //maxshift2[1]=1;
   //maxshift2[2]=1;
   //maxrot2[0]=1;
   //maxrot2[1]=1;
   //maxrot2[2]=1;

   TH2F *hpxvsiter =  new TH2F("hpxvsiter","Iteration vs Px;Iteration;Px",
     nitermost,0,nitermost,nbin2iter, -maxshift[0], maxshift[0]);
   hstore->Add(hpxvsiter);
   TH2F *hpyvsiter =  new TH2F("hpyvsiter","Iteration vs Py;Iteration;Py",
     nitermost,0,nitermost,nbin2iter, -maxshift[1], maxshift[1]);
   hstore->Add(hpyvsiter);
   TH2F *hpzvsiter =  new TH2F("hpzvsiter","Iteration vs Pz;Iteration;Pz",
     nitermost,0,nitermost,nbin2iter, -maxshift[2], maxshift[2]);
   hstore->Add(hpzvsiter);
   TH2F *hpavsiter =  new TH2F("hpavsiter","Iteration vs Pa;Iteration;Pa",
     nitermost,0,nitermost,nbin2iter, -maxrot[0], maxrot[0]);
   hstore->Add(hpavsiter);
   TH2F *hpbvsiter =  new TH2F("hpbvsiter","Iteration vs Pb;Iteration;Pb",
     nitermost,0,nitermost,nbin2iter, -maxrot[1], maxrot[1]);
   hstore->Add(hpbvsiter);
   TH2F *hpcvsiter =  new TH2F("hpcvsiter","Iteration vs Pc;Iteration;Pc",
     nitermost,0,nitermost,nbin2iter, -maxrot[2], maxrot[2]);
   hstore->Add(hpcvsiter);


   hpxvsiter->SetStats(kFALSE);
   hpyvsiter->SetStats(kFALSE);
   hpzvsiter->SetStats(kFALSE);
   hpavsiter->SetStats(kFALSE);
   hpbvsiter->SetStats(kFALSE);
   hpcvsiter->SetStats(kFALSE);

   can->cd(1); hpxvsiter->Fill(0.,0.,1.); hpxvsiter->Draw("");
   can->cd(3); hpxvsiter->Fill(0.,0.,1.); hpyvsiter->Draw("");
   can->cd(5); hpxvsiter->Fill(0.,0.,1.); hpzvsiter->Draw("");
   can->cd(2); hpxvsiter->Fill(0.,0.,1.); hpavsiter->Draw("");
   can->cd(4); hpxvsiter->Fill(0.,0.,1.); hpbvsiter->Draw("");
   can->cd(6); hpxvsiter->Fill(0.,0.,1.); hpcvsiter->Draw("");

   for (Int_t i=0;i<ALIMAX;i++) { 
     can->cd(1);
     TPolyLine* plx = (TPolyLine*) (thePolylinePxArr->At(i));
     if (plx!=0) { plx->SetLineWidth(1); plx->Draw();}
     can->cd(3);
     TPolyLine* ply = (TPolyLine*) (thePolylinePyArr->At(i));
     if (ply!=0) { ply->SetLineWidth(1); ply->Draw();}
     can->cd(5);
     TPolyLine* plz = (TPolyLine*) (thePolylinePzArr->At(i));
     if (plz!=0) { plz->SetLineWidth(1); plz->Draw();}
     can->cd(2);
     TPolyLine* pla = (TPolyLine*) (thePolylinePaArr->At(i));
     if (pla!=0) { pla->SetLineWidth(1); pla->Draw();}
     can->cd(4);
     TPolyLine* plb = (TPolyLine*) (thePolylinePbArr->At(i));
     if (plb!=0) { plb->SetLineWidth(1); plb->Draw();}
     can->cd(6);
     TPolyLine* plc = (TPolyLine*) (thePolylinePcArr->At(i));
     if (plc!=0) { plc->SetLineWidth(1); plc->Draw();}
   }

   can->Update();
   can->Print(psname,"ps"); 

   // -------------------------------------------------------------------------
   // shifts vs position

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); PlotTLHisto(hdxvslay,"BOX");
   can->cd(2); PlotTLHisto(hdyvslay,"BOX");
   can->cd(3); PlotTLHisto(hdzvslay,"BOX");
   can->Update();
   can->Print(psname,"ps"); 

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdxvsz->Draw("BOX");
   can->cd(2); hdyvsz->Draw("BOX");
   can->cd(3); hdzvsz->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdxvsr->Draw("BOX");
   can->cd(2); hdyvsr->Draw("BOX");
   can->cd(3); hdzvsr->Draw("BOX");
   //
   can->cd(4); hdphivsr->Draw("BOX");
   //
   can->Update();
   can->Print(psname,"ps"); 

   //
   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdphivsrl1->Draw("BOX");
   can->cd(2); hdphivsrl2->Draw("BOX");
   can->cd(3); hdphivsrl3->Draw("BOX");
   can->cd(4); hdphivsrl4->Draw("BOX");
   can->Update();
   can->Print(psname,"ps");
   //

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdxvsph->Draw("BOX");
   can->cd(2); hdyvsph->Draw("BOX");
   can->cd(3); hdzvsph->Draw("BOX");
   can->cd(4); hdxlocvsph->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 

   hdxvshit->SetStats(kFALSE);
   hdyvshit->SetStats(kFALSE);
   hdzvshit->SetStats(kFALSE);

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdxvshit->Draw("BOX");
   can->cd(2); hdyvshit->Draw("BOX");
   can->cd(3); hdzvshit->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 


   // -------------------------------------------------------------------------
   // rotations vs position

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); PlotTLHisto(hdavslay,"BOX");
   can->cd(2); PlotTLHisto(hdbvslay,"BOX");
   can->cd(3); PlotTLHisto(hdcvslay,"BOX");
   can->Update();
   can->Print(psname,"ps"); 

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdavsz->Draw("BOX");
   can->cd(2); hdbvsz->Draw("BOX");
   can->cd(3); hdcvsz->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdavsr->Draw("BOX");
   can->cd(2); hdbvsr->Draw("BOX");
   can->cd(3); hdcvsr->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,2);
   can->cd(1); hdavsph->Draw("BOX");
   can->cd(2); hdbvsph->Draw("BOX");
   can->cd(3); hdcvsph->Draw("BOX");
   can->Update();
   can->Print(psname,"ps"); 

   // parameter residuals -----------------------------------------------------

   can->Clear();
   globalTitle(isempty);
   can->Divide(2,3);

   can->cd(1); hresp1->Draw("");
   can->cd(3); hresp2->Draw("");
   can->cd(5); hresp3->Draw("");
   can->cd(2); hresp4->Draw("");
   can->cd(4); hresp5->Draw("");
   can->cd(6); hresp6->Draw("");

   can->Update();
   can->Print(psname,"ps"); 

   // correlations ----------------------------------------------------------
   can->Clear();
   globalTitle(isempty);
   can->Divide(3,3);
   
   int ipad = 1;
   // Check if we need to plot correlations and draw
   can->cd(ipad);
   if (floatVar.Contains("x") && floatVar.Contains("y")) {
     hdxvsdy->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("x") && floatVar.Contains("z")) {
     hdxvsdz->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("x") && floatVar.Contains("a")) {
     hdxvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("x") && floatVar.Contains("b")) {
     hdxvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("x") && floatVar.Contains("g")) {
     hdxvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("z")) {
     hdyvsdz->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("y") && floatVar.Contains("a")) {
     hdyvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("b")) {
     hdyvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("g")) {
     hdyvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("a")) {
     hdzvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("b")) {
     hdzvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("g")) {
     hdzvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("a") && floatVar.Contains("b")) {
     hdavsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("a") && floatVar.Contains("g")) {
     hdavsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("b") && floatVar.Contains("g")) {
     hdbvsdc->Draw("BOX");   can->cd(++ipad);
   } 

   can->Update();
   can->Print(psname,"ps");
  
   can->Clear();
   globalTitle(isempty);
   can->Divide(3,3);
   ipad = 1;
   // Check if we need to plot correlations and draw
   can->cd(ipad);
   if (floatVar.Contains("x") && floatVar.Contains("y")) {
     hstepdxvsdy->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("x") && floatVar.Contains("z")) {
     hstepdxvsdz->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("x") && floatVar.Contains("a")) {
     hstepdxvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("x") && floatVar.Contains("b")) {
     hstepdxvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("x") && floatVar.Contains("g")) {
     hstepdxvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("z")) {
     hstepdyvsdz->Draw("BOX");   can->cd(++ipad);
   } 
   if (floatVar.Contains("y") && floatVar.Contains("a")) {
     hstepdyvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("b")) {
     hstepdyvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("y") && floatVar.Contains("g")) {
     hstepdyvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("a")) {
     hstepdzvsda->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("b")) {
     hstepdzvsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("z") && floatVar.Contains("g")) {
     hstepdzvsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("a") && floatVar.Contains("b")) {
     hstepdavsdb->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("a") && floatVar.Contains("g")) {
     hstepdavsdc->Draw("BOX");   can->cd(++ipad);
   }
   if (floatVar.Contains("b") && floatVar.Contains("g")) {
     hstepdbvsdc->Draw("BOX");   can->cd(++ipad);
   } 
 
   can->Update();
   can->Print(psname,"ps");
}

//-----------------------------------------------------------------------------

void AnalysisLayers::globalTitle(TString s)
{
  if (s!="") {
    TPaveLabel* label = new TPaveLabel(0.3,0.99,0.7,1.03,s,"NDC");
    label->SetBorderSize(0);
    label->Draw();
  }
}

#endif
