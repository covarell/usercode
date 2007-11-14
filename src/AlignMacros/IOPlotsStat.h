
#ifndef IOPlotsStat_h
#define IOPlotsStat_h

#include "IOCombined.h"

class IOPlotsStat : public IOCombined {

 public:

  // constructor
  IOPlotsStat() {};

  // destructor
  ~IOPlotsStat() {};

  // book histos
  void analysis(void);

 private:

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void IOPlotsStat::analysis(void)
{

  // book ---------------------------------------------------------------------

  Float_t zmax=300., rmax=150., phmax=3.14159;
  Int_t nbin2=8;

  TH1F* alldets = new TH1F("alldets","Test",54,-27,27);
  TH1F* hitdets = new TH1F("hitdets","Test",54,-27,27);
  TH1F* ratio =  new TH1F("ratio","Ratio",54,-27,27);

  TH2F* alldetsvsz = new TH2F("alldetsvsz","Test",54,-27,27,nbin2,-zmax,zmax);
  TH2F* hitdetsvsz = new TH2F("hitdetsvsz","Test",54,-27,27,nbin2,-zmax,zmax);
  TH2F* ratioz     = new TH2F("ratioz","Test",54,-27,27,nbin2,-zmax,zmax);

  TH2F* alldetsvsr = new TH2F("alldetsvsr","Test",54,-27,27,nbin2,0,rmax);
  TH2F* hitdetsvsr = new TH2F("hitdetsvsr","Test",54,-27,27,nbin2,0,rmax);
  TH2F* ratior     = new TH2F("ratior","Test",54,-27,27,nbin2,0,rmax);

  TH2F* alldetsvsph = new TH2F("alldetsvsph","Test",54,-27,27,nbin2,-phmax,phmax);
  TH2F* hitdetsvsph = new TH2F("hitdetsvsph","Test",54,-27,27,nbin2,-phmax,phmax);
  TH2F* ratioph     = new TH2F("ratioph","Test",54,-27,27,nbin2,-phmax,phmax);



  static const int nnmin=5;

  //TString ele="Det";
  //static const int nmin[5]={5,10,15,30,50};   // dets
  //TString ele="Rod";
  //static const int nmin[5]={10,20,50,200,400};   // rods
  TString ele="Layer";
  static const int nmin[5]={500,750,1000,1500,2000};   // layers


  TString tit="Fraction of "+ele+"'s with N_{hit}>N_{min}";
  TString atit="Min. Hits per "+ele+" N_{min}";

  Float_t allarr[54][nnmin];
  Float_t hitarr[54][nnmin];
  Float_t ratarr[54][nnmin];

  for (Int_t i=0;i<54;i++) for(Int_t j=0;j<nnmin;j++) {
    hitarr[i][j]=0.; allarr[i][j]=0.;
  }

  TH2F* ratarrh = new TH2F("ratarrh",tit,54,-27,27,nnmin,0,nnmin);

  for (Int_t i=0;i<nnmin;i++) {
    Char_t slabel[4]; sprintf(slabel,"%d",nmin[i]);
    ratarrh->GetYaxis()->SetBinLabel(i+1,slabel);
  }
  ratarrh->GetYaxis()->SetTitle(atit);
  ratarrh->GetXaxis()->SetTitle("Layer");



  // table per components

  TString avtit="Average hit multiplicity per "+ele;
  TH1F* avdets = new TH1F("avdets",avtit,54,-27,27);
  avdets->GetYaxis()->SetTitle(avtit);
  avdets->GetXaxis()->SetTitle("Layer");

  Float_t ncomp[54]; // number of elements per component
  Float_t hitcomp[54]; // sum of hits
  Float_t avhits[54]; // average hits per component

  for  (Int_t i=0; i<54; i++) {
    ncomp[i]=0;
    hitcomp[i]=0;
    avhits[i]=0;
  }

  // fill ---------------------------------------------------------------------

  for (Int_t ia=0;ia<NAli;ia++) {
    if (AliUse[ia]==1) {

      Int_t tl=getlayer(ia);
      Float_t zpos=AliPosZ[ia];
      Float_t rpos=AliPosR[ia];
      Float_t phpos=AliPosPhi[ia];

      Int_t ci = tl+27;
      //Int_t ci=TMath::Abs(AliType[ia]);
      ncomp[ci]++;
      hitcomp[ci]+=AliNhit[ia];



      for (Int_t inm=0;inm<nnmin;inm++) { 
        allarr[tl+27][inm]++;
        if (AliNhit[ia]>nmin[inm]) hitarr[tl+27][inm]++;
      }

      alldets->Fill(tl,1.);
      alldetsvsz->Fill(tl,zpos,1.);
      alldetsvsr->Fill(tl,rpos,1.);
      alldetsvsph->Fill(tl,phpos,1.);
      if (AliNhit[ia]>10)  { 
	hitdets->Fill(tl,1.);
	hitdetsvsz->Fill(tl,zpos,1.);
	hitdetsvsr->Fill(tl,rpos,1.);
	hitdetsvsph->Fill(tl,phpos,1.);
      }
    }
  }

  ratio->Divide(hitdets,alldets);
  ratior->Divide(hitdetsvsr,alldetsvsr);
  ratioz->Divide(hitdetsvsz,alldetsvsz);
  ratioph->Divide(hitdetsvsph,alldetsvsph);

  for (Int_t il=0; il<54; il++) {
    for (Int_t inm=0;inm<nnmin;inm++) { 
      ratarr[il][inm]=0.;
      if (allarr[il][inm]>0.) ratarr[il][inm]=hitarr[il][inm]/allarr[il][inm];
      // printf("layer, ihit, content: %d %d %f\n",il,inm,ratarr[il][inm]);
      ratarrh->Fill(il-27,inm,ratarr[il][inm]);
    }
  }


  for (Int_t i=0;i<54;i++) {
    avhits[i]=hitcomp[i]/ncomp[i];
    avdets->Fill(i-27,avhits[i]);
    printf("layer, av hits %d %f \n",i,avhits[i]);
  }

  // plot ---------------------------------------------------------------------



   gSystem->Exec("rm -f IOPlotsStat.ps");
   TCanvas *can = new TCanvas("can", "Test");
   can->Print("IOPlotsStat.ps[","ps"); 

   // -------------------------------------------------------------------------

   can->Clear();
   can->Divide(1,1);

   can->cd(1); PlotTLHisto(ratarrh,"COLZ");

   can->Update();
   can->Print("IOPlotsStat.ps","ps"); 

   can->Clear();
   can->Divide(1,1);

   can->cd(1); PlotTLHisto(avdets);

   //can->cd(1); PlotTLHisto(ratioz,"COLZ");
   //can->cd(2); PlotTLHisto(ratior,"COLZ");
   //   can->cd(3); PlotTLHisto(ratioph,"COLZ");

   //   can->cd(2); PlotTLHisto(ratio);
   //can->cd(3); PlotTLHisto(ratio);

   can->Update();
   can->Print("IOPlotsStat.ps","ps"); 

   // -------------------------------------------------------------------------

   can->Print("IOPlotsStat.ps]","ps"); 
}



#endif
