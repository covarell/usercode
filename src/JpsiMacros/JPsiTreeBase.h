//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 22 16:30:41 2009 by ROOT version 5.18/00a
// from TTree T1/CMSSW Quarkonia tree
// found on file: /data/covarell/BtoJpsi_pt25/inclBtoJpsiMuMu_10TeV_analysis_2210_84.root
//////////////////////////////////////////////////////////

#ifndef JPsiTreeBase_h
#define JPsiTreeBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>

static const unsigned int MAXMC = 10;
static const unsigned int MAXMCMU = 100;
static const unsigned int MAXMU = 40;
static const unsigned int MAXMUCALO = 200;
static const unsigned int MAXQQ = 3000;
static const unsigned int MAXCHIC = 3000;
static const unsigned int MAXPRIVTX = 30;
static const unsigned int MAXTRIG = 10;

static const bool isAOD = true;
static const bool isOctX = true;

class JPsiTreeBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Mc_ProcessId;
   Double_t        Mc_EventScale;
   Double_t        Mc_EventWeight;
   Int_t           Mc_QQ_size;
   TClonesArray    *Mc_QQ_4mom;
   TClonesArray    *Mc_QQ_3vec;
   TClonesArray    *Mc_QQmoth_4mom;
   TClonesArray    *Mc_QQmoth_3vec;
   Int_t           Mc_QQmoth_id[MAXMC];   //[Mc_QQ_size]
   Int_t           Mc_QQmupl_indx[MAXMC];   //[Mc_QQ_size]
   Int_t           Mc_QQmumi_indx[MAXMC];   //[Mc_QQ_size]
   Int_t           Mc_mu_size;
   TClonesArray    *Mc_mu_4mom;
   TClonesArray    *Mc_mu_3vec;
   Int_t           Mc_mu_id[MAXMCMU];   //[Mc_mu_size]
   Int_t           Mc_mumoth_id[MAXMCMU];   //[Mc_mu_size]
   Int_t           Reco_mu_glb_size;
   TClonesArray    *Reco_mu_glb_4mom;
   TClonesArray    *Reco_mu_glb_track4mom;
   TClonesArray    *Reco_mu_glb_3vec;
   Double_t        Reco_mu_glb_phiErr[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_etaErr[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_ptErr[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_d0[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_d0err[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_dz[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_dzerr[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_normChi2[MAXMU];   //[Reco_mu_glb_size]
   Int_t           Reco_mu_glb_nhitstrack[MAXMU];   //[Reco_mu_glb_size]
   Int_t           Reco_mu_glb_nhitsDT[MAXMU];   //[Reco_mu_glb_size]
   Int_t           Reco_mu_glb_nhitsCSC[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_caloComp[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_segmComp[MAXMU];   //[Reco_mu_glb_size]
   Double_t        Reco_mu_glb_iso[MAXMU];   //[Reco_mu_glb_size]
   Int_t           Reco_mu_glb_charge[MAXMU];   //[Reco_mu_glb_size]
   Int_t           Reco_mu_trk_size;
   TClonesArray    *Reco_mu_trk_4mom;
   TClonesArray    *Reco_mu_trk_3vec;
   Double_t        Reco_mu_trk_phiErr[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_etaErr[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_ptErr[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_d0[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_d0err[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_dz[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_dzerr[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_normChi2[MAXMU];   //[Reco_mu_trk_size]
   Int_t           Reco_mu_trk_nhitstrack[MAXMU];   //[Reco_mu_trk_size]
   // Int_t           Reco_mu_trk_nhitsDT[MAXMU];   //[Reco_mu_trk_size]
   // Int_t           Reco_mu_trk_nhitsCSC[MAXMU];   //[Reco_mu_trk_size]
   Int_t           Reco_mu_trk_PIDmask[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_caloComp[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_segmComp[MAXMU];   //[Reco_mu_trk_size]
   Double_t        Reco_mu_trk_iso[MAXMU];   //[Reco_mu_trk_size]
   Int_t           Reco_mu_trk_charge[MAXMU];   //[Reco_mu_trk_size]
   Int_t           Reco_mu_cal_size;
   TClonesArray    *Reco_mu_cal_4mom;
   TClonesArray    *Reco_mu_cal_3vec;
   Double_t        Reco_mu_cal_phiErr[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_etaErr[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_ptErr[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_d0[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_d0err[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_dz[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_dzerr[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_normChi2[MAXMUCALO];   //[Reco_mu_cal_size]
   Int_t           Reco_mu_cal_nhitstrack[MAXMUCALO];   //[Reco_mu_cal_size]
   Double_t        Reco_mu_cal_caloComp[MAXMUCALO];   //[Reco_mu_cal_size]
   Int_t           Reco_mu_cal_charge[MAXMUCALO];   //[Reco_mu_cal_size]
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[MAXQQ];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   Int_t           Reco_QQ_mupl[MAXQQ];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi[MAXQQ];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mulpt[MAXQQ];   //[Reco_QQ_size]
   Int_t           Reco_QQ_muhpt[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_DeltaR[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_cosTheta[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_s[MAXQQ];   //[Reco_QQ_size]
   Char_t          Reco_QQ_VtxIsVal[MAXQQ];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_Vtx;
   Double_t        Reco_QQ_VxxE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_VyyE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_VzzE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_VyxE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_VzxE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_VzyE[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_lxy[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_lxyErr[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_normChi2[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_cosAlpha[MAXQQ];   //[Reco_QQ_size]
   Double_t        Reco_QQ_ctau[MAXQQ];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[MAXQQ];   //[Reco_QQ_size]
   Int_t           Reco_Chic_size;
   TClonesArray    *Reco_Chic_4mom;
   Int_t           Reco_Chic_OniaDaug[MAXCHIC];   //[Reco_Chic_size]
   Int_t           Reco_Chic_GammaDaug[MAXCHIC];   //[Reco_Chic_size]
   Double_t        Reco_Chic_DeltaM[MAXCHIC];   //[Reco_Chic_size]
   Double_t        Reco_BeamSpot_x;
   Double_t        Reco_BeamSpot_y;
   Double_t        Reco_BeamSpot_z;
   Double_t        Reco_BeamSpot_xxE;
   Double_t        Reco_BeamSpot_yyE;
   Double_t        Reco_BeamSpot_zzE;
   Double_t        Reco_BeamSpot_yxE;
   Double_t        Reco_BeamSpot_zyE;
   Double_t        Reco_BeamSpot_zxE;
   Int_t           Reco_PriVtx_size;
   TClonesArray    *Reco_PriVtx_3vec;
   Double_t        Reco_PriVtx_xxE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_yyE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_zzE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_yxE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_zyE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_zxE[MAXPRIVTX];   //[Reco_PriVtx_size]
   Int_t           Reco_PriVtx_trkSize[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_chi2[MAXPRIVTX];   //[Reco_PriVtx_size]
   Double_t        Reco_PriVtx_ndof[MAXPRIVTX];   //[Reco_PriVtx_size]
   Int_t           L1TBits_size;
   Char_t          L1TBits_accept[128];   //[L1TBits_size]
   Char_t          L1TGlobal_Decision;
   Int_t           L1_mu_size;
   TClonesArray    *L1_mu_4mom;
   Int_t           L1_mu_charge[50];   //[L1_mu_size]
   Int_t           HLTBits_size;
   Char_t          HLTBits_wasrun[8];   //[HLTBits_size]
   Char_t          HLTBits_accept[8];   //[HLTBits_size]
   Char_t          HLTBits_error[8];   //[HLTBits_size]
   Char_t          HLTGlobal_wasrun;
   Char_t          HLTGlobal_Decision;
   Char_t          HLTGlobal_error;
   Int_t           HLT1Mu3_L3_size;
   TClonesArray    *HLT1Mu3_L3_4mom;
   Int_t           HLT1Mu3_L3_id[MAXTRIG];   //[HLT1Mu3_L3_size]
   Int_t           HLT1Mu5_L3_size;
   TClonesArray    *HLT1Mu5_L3_4mom;
   Int_t           HLT1Mu5_L3_id[MAXTRIG];   //[HLT1Mu5_L3_size]
   Int_t           HLT1Mu9_L3_size;
   TClonesArray    *HLT1Mu9_L3_4mom;
   Int_t           HLT1Mu9_L3_id[MAXTRIG];   //[HLT1Mu9_L3_size]
   Int_t           HLT1Mu11_L3_size;
   TClonesArray    *HLT1Mu11_L3_4mom;
   Int_t           HLT1Mu11_L3_id[MAXTRIG];   //[HLT1Mu11_L3_size]
   Int_t           HLT2IsoMu3_L3_size;
   TClonesArray    *HLT2IsoMu3_L3_4mom;
   Int_t           HLT2IsoMu3_L3_id[MAXTRIG];   //[HLT2IsoMu3_L3_size]
   Int_t           HLT2Mu3_L3_size;
   TClonesArray    *HLT2Mu3_L3_4mom;
   Int_t           HLT2Mu3_L3_id[MAXTRIG];   //[HLT2Mu3_L3_size]
   Int_t           HLTJpsi2Mu_L3_size;
   TClonesArray    *HLTJpsi2Mu_L3_4mom;
   Int_t           HLTJpsi2Mu_L3_id[MAXTRIG];   //[HLTJpsi2Mu_L3_size]
   Int_t           HLTUpsilon2Mu_L3_size;
   TClonesArray    *HLTUpsilon2Mu_L3_4mom;
   Int_t           HLTUpsilon2Mu_L3_id[MAXTRIG];   //[HLTUpsilon2Mu_L3_size]

   // List of branches
   TBranch        *b_Mc_ProcessId;   //!
   TBranch        *b_Mc_EventScale;   //!
   TBranch        *b_Mc_EventWeight;   //!
   TBranch        *b_Mc_QQ_size;   //!
   TBranch        *b_Mc_QQ_4mom;   //!
   TBranch        *b_Mc_QQ_3vec;   //!
   TBranch        *b_Mc_QQmoth_4mom;   //!
   TBranch        *b_Mc_QQmoth_3vec;   //!
   TBranch        *b_Mc_QQmoth_id;   //!
   TBranch        *b_Mc_QQmupl_indx;   //!
   TBranch        *b_Mc_QQmumi_indx;   //!
   TBranch        *b_Mc_mu_size;   //!
   TBranch        *b_Mc_mu_4mom;   //!
   TBranch        *b_Mc_mu_3vec;   //!
   TBranch        *b_Mc_mu_id;   //!
   TBranch        *b_Mc_mumoth_id;   //!
   TBranch        *b_Reco_mu_glb_size;   //!
   TBranch        *b_Reco_mu_glb_4mom;   //!
   TBranch        *b_Reco_mu_glb_track4mom;   //!
   TBranch        *b_Reco_mu_glb_3vec;   //!
   TBranch        *b_Reco_mu_glb_phiErr;   //!
   TBranch        *b_Reco_mu_glb_etaErr;   //!
   TBranch        *b_Reco_mu_glb_ptErr;   //!
   TBranch        *b_Reco_mu_glb_d0;   //!
   TBranch        *b_Reco_mu_glb_d0err;   //!
   TBranch        *b_Reco_mu_glb_dz;   //!
   TBranch        *b_Reco_mu_glb_dzerr;   //!
   TBranch        *b_Reco_mu_glb_normChi2;   //!
   TBranch        *b_Reco_mu_glb_nhitstrack;   //!
   TBranch        *b_Reco_mu_glb_nhitsDT;   //!
   TBranch        *b_Reco_mu_glb_nhitsCSC;   //!
   TBranch        *b_Reco_mu_glb_caloComp;   //!
   TBranch        *b_Reco_mu_glb_segmComp;   //!
   TBranch        *b_Reco_mu_glb_iso;   //!
   TBranch        *b_Reco_mu_glb_charge;   //!
   TBranch        *b_Reco_mu_trk_size;   //!
   TBranch        *b_Reco_mu_trk_4mom;   //!
   TBranch        *b_Reco_mu_trk_3vec;   //!
   TBranch        *b_Reco_mu_trk_phiErr;   //!
   TBranch        *b_Reco_mu_trk_etaErr;   //!
   TBranch        *b_Reco_mu_trk_ptErr;   //!
   TBranch        *b_Reco_mu_trk_d0;   //!
   TBranch        *b_Reco_mu_trk_d0err;   //!
   TBranch        *b_Reco_mu_trk_dz;   //!
   TBranch        *b_Reco_mu_trk_dzerr;   //!
   TBranch        *b_Reco_mu_trk_normChi2;   //!
   TBranch        *b_Reco_mu_trk_nhitstrack;   //!
   // TBranch        *b_Reco_mu_trk_nhitsDT;   //!
   // TBranch        *b_Reco_mu_trk_nhitsCSC;   //!
   TBranch        *b_Reco_mu_trk_PIDmask;   //!
   TBranch        *b_Reco_mu_trk_caloComp;   //!
   TBranch        *b_Reco_mu_trk_segmComp;   //!
   TBranch        *b_Reco_mu_trk_iso;   //!
   TBranch        *b_Reco_mu_trk_charge;   //!
   TBranch        *b_Reco_mu_cal_size;   //!
   TBranch        *b_Reco_mu_cal_4mom;   //!
   TBranch        *b_Reco_mu_cal_3vec;   //!
   TBranch        *b_Reco_mu_cal_phiErr;   //!
   TBranch        *b_Reco_mu_cal_etaErr;   //!
   TBranch        *b_Reco_mu_cal_ptErr;   //!
   TBranch        *b_Reco_mu_cal_d0;   //!
   TBranch        *b_Reco_mu_cal_d0err;   //!
   TBranch        *b_Reco_mu_cal_dz;   //!
   TBranch        *b_Reco_mu_cal_dzerr;   //!
   TBranch        *b_Reco_mu_cal_normChi2;   //!
   TBranch        *b_Reco_mu_cal_nhitstrack;   //!
   TBranch        *b_Reco_mu_cal_caloComp;   //!
   TBranch        *b_Reco_mu_cal_charge;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl;   //!
   TBranch        *b_Reco_QQ_mumi;   //!
   TBranch        *b_Reco_QQ_mulpt;   //!
   TBranch        *b_Reco_QQ_muhpt;   //!
   TBranch        *b_Reco_QQ_DeltaR;   //!
   TBranch        *b_Reco_QQ_cosTheta;   //!
   TBranch        *b_Reco_QQ_s;   //!
   TBranch        *b_Reco_QQ_VtxIsVal;   //!
   TBranch        *b_Reco_QQ_Vtx;   //!
   TBranch        *b_Reco_QQ_VxxE;   //!
   TBranch        *b_Reco_QQ_VyyE;   //!
   TBranch        *b_Reco_QQ_VzzE;   //!
   TBranch        *b_Reco_QQ_VyxE;   //!
   TBranch        *b_Reco_QQ_VzxE;   //!
   TBranch        *b_Reco_QQ_VzyE;   //!
   TBranch        *b_Reco_QQ_lxy;   //!
   TBranch        *b_Reco_QQ_lxyErr;   //!
   TBranch        *b_Reco_QQ_normChi2;   //!
   TBranch        *b_Reco_QQ_cosAlpha;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_Chic_size;   //!
   TBranch        *b_Reco_Chic_4mom;   //!
   TBranch        *b_Reco_Chic_OniaDaug;   //!
   TBranch        *b_Reco_Chic_GammaDaug;   //!
   TBranch        *b_Reco_Chic_DeltaM;   //!
   TBranch        *b_Reco_BeamSpot_x;   //!
   TBranch        *b_Reco_BeamSpot_y;   //!
   TBranch        *b_Reco_BeamSpot_z;   //!
   TBranch        *b_Reco_BeamSpot_xxE;   //!
   TBranch        *b_Reco_BeamSpot_yyE;   //!
   TBranch        *b_Reco_BeamSpot_zzE;   //!
   TBranch        *b_Reco_BeamSpot_yxE;   //!
   TBranch        *b_Reco_BeamSpot_zyE;   //!
   TBranch        *b_Reco_BeamSpot_zxE;   //!
   TBranch        *b_Reco_PriVtx_size;   //!
   TBranch        *b_Reco_PriVtx_3vec;   //!
   TBranch        *b_Reco_PriVtx_xxE;   //!
   TBranch        *b_Reco_PriVtx_yyE;   //!
   TBranch        *b_Reco_PriVtx_zzE;   //!
   TBranch        *b_Reco_PriVtx_yxE;   //!
   TBranch        *b_Reco_PriVtx_zyE;   //!
   TBranch        *b_Reco_PriVtx_zxE;   //!
   TBranch        *b_Reco_PriVtx_trkSize;   //!
   TBranch        *b_Reco_PriVtx_chi2;   //!
   TBranch        *b_Reco_PriVtx_ndof;   //!
   TBranch        *b_L1TBits_size;   //!
   TBranch        *b_L1TBits_accept;   //!
   TBranch        *b_L1TGlobal_Decision;   //!
   TBranch        *b_L1_mu_size;   //!
   TBranch        *b_L1_mu_4mom;   //!
   TBranch        *b_L1_mu_charge;   //!
   TBranch        *b_HLTBits_size;   //!
   TBranch        *b_HLTBits_wasrun;   //!
   TBranch        *b_HLTBits_accept;   //!
   TBranch        *b_HLTBits_error;   //!
   TBranch        *b_HLTGlobal_wasrun;   //!
   TBranch        *b_HLTGlobal_Decision;   //!
   TBranch        *b_HLTGlobal_error;   //!
   TBranch        *b_HLT1Mu3_L3_size;   //!
   TBranch        *b_HLT1Mu3_L3_4mom;   //!
   TBranch        *b_HLT1Mu3_L3_id;   //!
   TBranch        *b_HLT1Mu5_L3_size;   //!
   TBranch        *b_HLT1Mu5_L3_4mom;   //!
   TBranch        *b_HLT1Mu5_L3_id;   //!
   TBranch        *b_HLT1Mu9_L3_size;   //!
   TBranch        *b_HLT1Mu9_L3_4mom;   //!
   TBranch        *b_HLT1Mu9_L3_id;   //!
   TBranch        *b_HLT1Mu11_L3_size;   //!
   TBranch        *b_HLT1Mu11_L3_4mom;   //!
   TBranch        *b_HLT1Mu11_L3_id;   //!
   TBranch        *b_HLT2IsoMu3_L3_size;   //!
   TBranch        *b_HLT2IsoMu3_L3_4mom;   //!
   TBranch        *b_HLT2IsoMu3_L3_id;   //!
   TBranch        *b_HLT2Mu3_L3_size;   //!
   TBranch        *b_HLT2Mu3_L3_4mom;   //!
   TBranch        *b_HLT2Mu3_L3_id;   //!
   TBranch        *b_HLTJpsi2Mu_L3_size;   //!
   TBranch        *b_HLTJpsi2Mu_L3_4mom;   //!
   TBranch        *b_HLTJpsi2Mu_L3_id;   //!
   TBranch        *b_HLTUpsilon2Mu_L3_size;   //!
   TBranch        *b_HLTUpsilon2Mu_L3_4mom;   //!
   TBranch        *b_HLTUpsilon2Mu_L3_id;   //!

   JPsiTreeBase(TTree *tree=0);
   virtual ~JPsiTreeBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JPsiTreeBase_cxx
JPsiTreeBase::JPsiTreeBase(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/covarell/BtoJpsi_pt25/inclBtoJpsiMuMu_10TeV_analysis_2210_84.root");
      if (!f) {
         f = new TFile("/data/covarell/BtoJpsi_pt25/inclBtoJpsiMuMu_10TeV_analysis_2210_84.root");
      }
      tree = (TTree*)gDirectory->Get("T1");

   }
   Init(tree);
}

JPsiTreeBase::~JPsiTreeBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JPsiTreeBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JPsiTreeBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JPsiTreeBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Mc_QQ_4mom = 0;
   Mc_QQ_3vec = 0;
   Mc_QQmoth_4mom = 0;
   Mc_QQmoth_3vec = 0;
   Mc_mu_4mom = 0;
   Mc_mu_3vec = 0;
   Reco_mu_glb_4mom = 0;
   Reco_mu_glb_track4mom = 0;
   Reco_mu_glb_3vec = 0;
   Reco_mu_trk_4mom = 0;
   Reco_mu_trk_3vec = 0;
   Reco_mu_cal_4mom = 0;
   Reco_mu_cal_3vec = 0;
   Reco_QQ_4mom = 0;
   Reco_QQ_Vtx = 0;
   Reco_Chic_4mom = 0;
   Reco_PriVtx_3vec = 0;
   L1_mu_4mom = 0;
   HLT1Mu3_L3_4mom = 0;
   HLT1Mu5_L3_4mom = 0;
   HLT1Mu9_L3_4mom = 0;
   HLT1Mu11_L3_4mom = 0;
   HLT2IsoMu3_L3_4mom = 0;
   HLT2Mu3_L3_4mom = 0;
   HLTJpsi2Mu_L3_4mom = 0;
   HLTUpsilon2Mu_L3_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Mc_ProcessId", &Mc_ProcessId, &b_Mc_ProcessId);
   fChain->SetBranchAddress("Mc_EventScale", &Mc_EventScale, &b_Mc_EventScale);
   fChain->SetBranchAddress("Mc_EventWeight", &Mc_EventWeight, &b_Mc_EventWeight);
   fChain->SetBranchAddress("Mc_QQ_size", &Mc_QQ_size, &b_Mc_QQ_size);
   fChain->SetBranchAddress("Mc_QQ_4mom", &Mc_QQ_4mom, &b_Mc_QQ_4mom);
   fChain->SetBranchAddress("Mc_QQ_3vec", &Mc_QQ_3vec, &b_Mc_QQ_3vec);
   fChain->SetBranchAddress("Mc_QQmoth_4mom", &Mc_QQmoth_4mom, &b_Mc_QQmoth_4mom);
   fChain->SetBranchAddress("Mc_QQmoth_3vec", &Mc_QQmoth_3vec, &b_Mc_QQmoth_3vec);
   fChain->SetBranchAddress("Mc_QQmoth_id", Mc_QQmoth_id, &b_Mc_QQmoth_id);
   fChain->SetBranchAddress("Mc_QQmupl_indx", Mc_QQmupl_indx, &b_Mc_QQmupl_indx);
   fChain->SetBranchAddress("Mc_QQmumi_indx", Mc_QQmumi_indx, &b_Mc_QQmumi_indx);
   fChain->SetBranchAddress("Mc_mu_size", &Mc_mu_size, &b_Mc_mu_size);
   fChain->SetBranchAddress("Mc_mu_4mom", &Mc_mu_4mom, &b_Mc_mu_4mom);
   fChain->SetBranchAddress("Mc_mu_3vec", &Mc_mu_3vec, &b_Mc_mu_3vec);
   fChain->SetBranchAddress("Mc_mu_id", Mc_mu_id, &b_Mc_mu_id);
   fChain->SetBranchAddress("Mc_mumoth_id", Mc_mumoth_id, &b_Mc_mumoth_id);
   fChain->SetBranchAddress("Reco_mu_glb_size", &Reco_mu_glb_size, &b_Reco_mu_glb_size);
   fChain->SetBranchAddress("Reco_mu_glb_4mom", &Reco_mu_glb_4mom, &b_Reco_mu_glb_4mom);
   fChain->SetBranchAddress("Reco_mu_glb_track4mom", &Reco_mu_glb_track4mom, &b_Reco_mu_glb_track4mom);
   fChain->SetBranchAddress("Reco_mu_glb_3vec", &Reco_mu_glb_3vec, &b_Reco_mu_glb_3vec);
   fChain->SetBranchAddress("Reco_mu_glb_phiErr", Reco_mu_glb_phiErr, &b_Reco_mu_glb_phiErr);
   fChain->SetBranchAddress("Reco_mu_glb_etaErr", Reco_mu_glb_etaErr, &b_Reco_mu_glb_etaErr);
   fChain->SetBranchAddress("Reco_mu_glb_ptErr", Reco_mu_glb_ptErr, &b_Reco_mu_glb_ptErr);
   fChain->SetBranchAddress("Reco_mu_glb_d0", Reco_mu_glb_d0, &b_Reco_mu_glb_d0);
   fChain->SetBranchAddress("Reco_mu_glb_d0err", Reco_mu_glb_d0err, &b_Reco_mu_glb_d0err);
   fChain->SetBranchAddress("Reco_mu_glb_dz", Reco_mu_glb_dz, &b_Reco_mu_glb_dz);
   fChain->SetBranchAddress("Reco_mu_glb_dzerr", Reco_mu_glb_dzerr, &b_Reco_mu_glb_dzerr);
   fChain->SetBranchAddress("Reco_mu_glb_normChi2", Reco_mu_glb_normChi2, &b_Reco_mu_glb_normChi2);
   fChain->SetBranchAddress("Reco_mu_glb_nhitstrack", Reco_mu_glb_nhitstrack, &b_Reco_mu_glb_nhitstrack);
   if (!isAOD) {
     fChain->SetBranchAddress("Reco_mu_glb_nhitsDT", Reco_mu_glb_nhitsDT, &b_Reco_mu_glb_nhitsDT);
     fChain->SetBranchAddress("Reco_mu_glb_nhitsCSC", Reco_mu_glb_nhitsCSC, &b_Reco_mu_glb_nhitsCSC);
   }
   fChain->SetBranchAddress("Reco_mu_glb_caloComp", Reco_mu_glb_caloComp, &b_Reco_mu_glb_caloComp);
   fChain->SetBranchAddress("Reco_mu_glb_segmComp", Reco_mu_glb_segmComp, &b_Reco_mu_glb_segmComp);
   fChain->SetBranchAddress("Reco_mu_glb_iso", Reco_mu_glb_iso, &b_Reco_mu_glb_iso);
   fChain->SetBranchAddress("Reco_mu_glb_charge", Reco_mu_glb_charge, &b_Reco_mu_glb_charge);
   fChain->SetBranchAddress("Reco_mu_trk_size", &Reco_mu_trk_size, &b_Reco_mu_trk_size);
   fChain->SetBranchAddress("Reco_mu_trk_4mom", &Reco_mu_trk_4mom, &b_Reco_mu_trk_4mom);
   fChain->SetBranchAddress("Reco_mu_trk_3vec", &Reco_mu_trk_3vec, &b_Reco_mu_trk_3vec);
   fChain->SetBranchAddress("Reco_mu_trk_phiErr", Reco_mu_trk_phiErr, &b_Reco_mu_trk_phiErr);
   fChain->SetBranchAddress("Reco_mu_trk_etaErr", Reco_mu_trk_etaErr, &b_Reco_mu_trk_etaErr);
   fChain->SetBranchAddress("Reco_mu_trk_ptErr", Reco_mu_trk_ptErr, &b_Reco_mu_trk_ptErr);
   fChain->SetBranchAddress("Reco_mu_trk_d0", Reco_mu_trk_d0, &b_Reco_mu_trk_d0);
   fChain->SetBranchAddress("Reco_mu_trk_d0err", Reco_mu_trk_d0err, &b_Reco_mu_trk_d0err);
   fChain->SetBranchAddress("Reco_mu_trk_dz", Reco_mu_trk_dz, &b_Reco_mu_trk_dz);
   fChain->SetBranchAddress("Reco_mu_trk_dzerr", Reco_mu_trk_dzerr, &b_Reco_mu_trk_dzerr);
   fChain->SetBranchAddress("Reco_mu_trk_normChi2", Reco_mu_trk_normChi2, &b_Reco_mu_trk_normChi2);
   fChain->SetBranchAddress("Reco_mu_trk_nhitstrack", Reco_mu_trk_nhitstrack, &b_Reco_mu_trk_nhitstrack);
   // fChain->SetBranchAddress("Reco_mu_trk_nhitsDT", Reco_mu_trk_nhitsDT, &b_Reco_mu_trk_nhitsDT);
   // fChain->SetBranchAddress("Reco_mu_trk_nhitsCSC", Reco_mu_trk_nhitsCSC, &b_Reco_mu_trk_nhitsCSC);
   fChain->SetBranchAddress("Reco_mu_trk_PIDmask", Reco_mu_trk_PIDmask, &b_Reco_mu_trk_PIDmask);
   fChain->SetBranchAddress("Reco_mu_trk_caloComp", Reco_mu_trk_caloComp, &b_Reco_mu_trk_caloComp);
   fChain->SetBranchAddress("Reco_mu_trk_segmComp", Reco_mu_trk_segmComp, &b_Reco_mu_trk_segmComp);
   fChain->SetBranchAddress("Reco_mu_trk_iso", Reco_mu_trk_iso, &b_Reco_mu_trk_iso);
   fChain->SetBranchAddress("Reco_mu_trk_charge", Reco_mu_trk_charge, &b_Reco_mu_trk_charge);
   if (!isOctX) {
     fChain->SetBranchAddress("Reco_mu_cal_size", &Reco_mu_cal_size, &b_Reco_mu_cal_size);
     fChain->SetBranchAddress("Reco_mu_cal_4mom", &Reco_mu_cal_4mom, &b_Reco_mu_cal_4mom);
     fChain->SetBranchAddress("Reco_mu_cal_3vec", &Reco_mu_cal_3vec, &b_Reco_mu_cal_3vec);
     fChain->SetBranchAddress("Reco_mu_cal_phiErr", Reco_mu_cal_phiErr, &b_Reco_mu_cal_phiErr);
     fChain->SetBranchAddress("Reco_mu_cal_etaErr", Reco_mu_cal_etaErr, &b_Reco_mu_cal_etaErr);
     fChain->SetBranchAddress("Reco_mu_cal_ptErr", Reco_mu_cal_ptErr, &b_Reco_mu_cal_ptErr);
     fChain->SetBranchAddress("Reco_mu_cal_d0", Reco_mu_cal_d0, &b_Reco_mu_cal_d0);
     fChain->SetBranchAddress("Reco_mu_cal_d0err", Reco_mu_cal_d0err, &b_Reco_mu_cal_d0err);
     fChain->SetBranchAddress("Reco_mu_cal_dz", Reco_mu_cal_dz, &b_Reco_mu_cal_dz);
     fChain->SetBranchAddress("Reco_mu_cal_dzerr", Reco_mu_cal_dzerr, &b_Reco_mu_cal_dzerr);
     fChain->SetBranchAddress("Reco_mu_cal_normChi2", Reco_mu_cal_normChi2, &b_Reco_mu_cal_normChi2);
     fChain->SetBranchAddress("Reco_mu_cal_nhitstrack", Reco_mu_cal_nhitstrack, &b_Reco_mu_cal_nhitstrack);
     fChain->SetBranchAddress("Reco_mu_cal_caloComp", Reco_mu_cal_caloComp, &b_Reco_mu_cal_caloComp);
     fChain->SetBranchAddress("Reco_mu_cal_charge", Reco_mu_cal_charge, &b_Reco_mu_cal_charge);
   }
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl", Reco_QQ_mupl, &b_Reco_QQ_mupl);
   fChain->SetBranchAddress("Reco_QQ_mumi", Reco_QQ_mumi, &b_Reco_QQ_mumi);
   fChain->SetBranchAddress("Reco_QQ_mulpt", Reco_QQ_mulpt, &b_Reco_QQ_mulpt);
   fChain->SetBranchAddress("Reco_QQ_muhpt", Reco_QQ_muhpt, &b_Reco_QQ_muhpt);
   fChain->SetBranchAddress("Reco_QQ_DeltaR", Reco_QQ_DeltaR, &b_Reco_QQ_DeltaR);
   fChain->SetBranchAddress("Reco_QQ_cosTheta", Reco_QQ_cosTheta, &b_Reco_QQ_cosTheta);
   fChain->SetBranchAddress("Reco_QQ_s", Reco_QQ_s, &b_Reco_QQ_s);
   fChain->SetBranchAddress("Reco_QQ_VtxIsVal", Reco_QQ_VtxIsVal, &b_Reco_QQ_VtxIsVal);
   fChain->SetBranchAddress("Reco_QQ_Vtx", &Reco_QQ_Vtx, &b_Reco_QQ_Vtx);
   fChain->SetBranchAddress("Reco_QQ_VxxE", Reco_QQ_VxxE, &b_Reco_QQ_VxxE);
   fChain->SetBranchAddress("Reco_QQ_VyyE", Reco_QQ_VyyE, &b_Reco_QQ_VyyE);
   fChain->SetBranchAddress("Reco_QQ_VzzE", Reco_QQ_VzzE, &b_Reco_QQ_VzzE);
   fChain->SetBranchAddress("Reco_QQ_VyxE", Reco_QQ_VyxE, &b_Reco_QQ_VyxE);
   fChain->SetBranchAddress("Reco_QQ_VzxE", Reco_QQ_VzxE, &b_Reco_QQ_VzxE);
   fChain->SetBranchAddress("Reco_QQ_VzyE", Reco_QQ_VzyE, &b_Reco_QQ_VzyE);
   fChain->SetBranchAddress("Reco_QQ_lxy", Reco_QQ_lxy, &b_Reco_QQ_lxy);
   fChain->SetBranchAddress("Reco_QQ_lxyErr", Reco_QQ_lxyErr, &b_Reco_QQ_lxyErr);
   fChain->SetBranchAddress("Reco_QQ_normChi2", Reco_QQ_normChi2, &b_Reco_QQ_normChi2);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   if (!isOctX) {
     fChain->SetBranchAddress("Reco_Chic_size", &Reco_Chic_size, &b_Reco_Chic_size);
     fChain->SetBranchAddress("Reco_Chic_4mom", &Reco_Chic_4mom, &b_Reco_Chic_4mom);
     fChain->SetBranchAddress("Reco_Chic_OniaDaug", Reco_Chic_OniaDaug, &b_Reco_Chic_OniaDaug);
     fChain->SetBranchAddress("Reco_Chic_GammaDaug", Reco_Chic_GammaDaug, &b_Reco_Chic_GammaDaug);
     fChain->SetBranchAddress("Reco_Chic_DeltaM", Reco_Chic_DeltaM, &b_Reco_Chic_DeltaM);
   }
   fChain->SetBranchAddress("Reco_BeamSpot_x", &Reco_BeamSpot_x, &b_Reco_BeamSpot_x);
   fChain->SetBranchAddress("Reco_BeamSpot_y", &Reco_BeamSpot_y, &b_Reco_BeamSpot_y);
   fChain->SetBranchAddress("Reco_BeamSpot_z", &Reco_BeamSpot_z, &b_Reco_BeamSpot_z);
   fChain->SetBranchAddress("Reco_BeamSpot_xxE", &Reco_BeamSpot_xxE, &b_Reco_BeamSpot_xxE);
   fChain->SetBranchAddress("Reco_BeamSpot_yyE", &Reco_BeamSpot_yyE, &b_Reco_BeamSpot_yyE);
   fChain->SetBranchAddress("Reco_BeamSpot_zzE", &Reco_BeamSpot_zzE, &b_Reco_BeamSpot_zzE);
   fChain->SetBranchAddress("Reco_BeamSpot_yxE", &Reco_BeamSpot_yxE, &b_Reco_BeamSpot_yxE);
   fChain->SetBranchAddress("Reco_BeamSpot_zyE", &Reco_BeamSpot_zyE, &b_Reco_BeamSpot_zyE);
   fChain->SetBranchAddress("Reco_BeamSpot_zxE", &Reco_BeamSpot_zxE, &b_Reco_BeamSpot_zxE);
   fChain->SetBranchAddress("Reco_PriVtx_size", &Reco_PriVtx_size, &b_Reco_PriVtx_size);
   fChain->SetBranchAddress("Reco_PriVtx_3vec", &Reco_PriVtx_3vec, &b_Reco_PriVtx_3vec);
   fChain->SetBranchAddress("Reco_PriVtx_xxE", Reco_PriVtx_xxE, &b_Reco_PriVtx_xxE);
   fChain->SetBranchAddress("Reco_PriVtx_yyE", Reco_PriVtx_yyE, &b_Reco_PriVtx_yyE);
   fChain->SetBranchAddress("Reco_PriVtx_zzE", Reco_PriVtx_zzE, &b_Reco_PriVtx_zzE);
   fChain->SetBranchAddress("Reco_PriVtx_yxE", Reco_PriVtx_yxE, &b_Reco_PriVtx_yxE);
   fChain->SetBranchAddress("Reco_PriVtx_zyE", Reco_PriVtx_zyE, &b_Reco_PriVtx_zyE);
   fChain->SetBranchAddress("Reco_PriVtx_zxE", Reco_PriVtx_zxE, &b_Reco_PriVtx_zxE);
   fChain->SetBranchAddress("Reco_PriVtx_trkSize", Reco_PriVtx_trkSize, &b_Reco_PriVtx_trkSize);
   fChain->SetBranchAddress("Reco_PriVtx_chi2", Reco_PriVtx_chi2, &b_Reco_PriVtx_chi2);
   fChain->SetBranchAddress("Reco_PriVtx_ndof", Reco_PriVtx_ndof, &b_Reco_PriVtx_ndof);
   fChain->SetBranchAddress("L1TBits_size", &L1TBits_size, &b_L1TBits_size);
   fChain->SetBranchAddress("L1TBits_accept", L1TBits_accept, &b_L1TBits_accept);
   fChain->SetBranchAddress("L1TGlobal_Decision", &L1TGlobal_Decision, &b_L1TGlobal_Decision);
   fChain->SetBranchAddress("L1_mu_size", &L1_mu_size, &b_L1_mu_size);
   fChain->SetBranchAddress("L1_mu_4mom", &L1_mu_4mom, &b_L1_mu_4mom);
   fChain->SetBranchAddress("L1_mu_charge", L1_mu_charge, &b_L1_mu_charge);
   fChain->SetBranchAddress("HLTBits_size", &HLTBits_size, &b_HLTBits_size);
   fChain->SetBranchAddress("HLTBits_wasrun", HLTBits_wasrun, &b_HLTBits_wasrun);
   fChain->SetBranchAddress("HLTBits_accept", HLTBits_accept, &b_HLTBits_accept);
   fChain->SetBranchAddress("HLTBits_error", HLTBits_error, &b_HLTBits_error);
   fChain->SetBranchAddress("HLTGlobal_wasrun", &HLTGlobal_wasrun, &b_HLTGlobal_wasrun);
   fChain->SetBranchAddress("HLTGlobal_Decision", &HLTGlobal_Decision, &b_HLTGlobal_Decision);
   fChain->SetBranchAddress("HLTGlobal_error", &HLTGlobal_error, &b_HLTGlobal_error);
   fChain->SetBranchAddress("HLT1Mu3_L3_size", &HLT1Mu3_L3_size, &b_HLT1Mu3_L3_size);
   fChain->SetBranchAddress("HLT1Mu3_L3_4mom", &HLT1Mu3_L3_4mom, &b_HLT1Mu3_L3_4mom);
   fChain->SetBranchAddress("HLT1Mu3_L3_id", HLT1Mu3_L3_id, &b_HLT1Mu3_L3_id);
   fChain->SetBranchAddress("HLT1Mu5_L3_size", &HLT1Mu5_L3_size, &b_HLT1Mu5_L3_size);
   fChain->SetBranchAddress("HLT1Mu5_L3_4mom", &HLT1Mu5_L3_4mom, &b_HLT1Mu5_L3_4mom);
   fChain->SetBranchAddress("HLT1Mu5_L3_id", HLT1Mu5_L3_id, &b_HLT1Mu5_L3_id);
   fChain->SetBranchAddress("HLT1Mu9_L3_size", &HLT1Mu9_L3_size, &b_HLT1Mu9_L3_size);
   fChain->SetBranchAddress("HLT1Mu9_L3_4mom", &HLT1Mu9_L3_4mom, &b_HLT1Mu9_L3_4mom);
   fChain->SetBranchAddress("HLT1Mu9_L3_id", HLT1Mu9_L3_id, &b_HLT1Mu9_L3_id);
   fChain->SetBranchAddress("HLT1Mu11_L3_size", &HLT1Mu11_L3_size, &b_HLT1Mu11_L3_size);
   fChain->SetBranchAddress("HLT1Mu11_L3_4mom", &HLT1Mu11_L3_4mom, &b_HLT1Mu11_L3_4mom);
   fChain->SetBranchAddress("HLT1Mu11_L3_id", HLT1Mu11_L3_id, &b_HLT1Mu11_L3_id);
   fChain->SetBranchAddress("HLT2IsoMu3_L3_size", &HLT2IsoMu3_L3_size, &b_HLT2IsoMu3_L3_size);
   fChain->SetBranchAddress("HLT2IsoMu3_L3_4mom", &HLT2IsoMu3_L3_4mom, &b_HLT2IsoMu3_L3_4mom);
   fChain->SetBranchAddress("HLT2IsoMu3_L3_id", HLT2IsoMu3_L3_id, &b_HLT2IsoMu3_L3_id);
   fChain->SetBranchAddress("HLT2Mu3_L3_size", &HLT2Mu3_L3_size, &b_HLT2Mu3_L3_size);
   fChain->SetBranchAddress("HLT2Mu3_L3_4mom", &HLT2Mu3_L3_4mom, &b_HLT2Mu3_L3_4mom);
   fChain->SetBranchAddress("HLT2Mu3_L3_id", HLT2Mu3_L3_id, &b_HLT2Mu3_L3_id);
   fChain->SetBranchAddress("HLTJpsi2Mu_L3_size", &HLTJpsi2Mu_L3_size, &b_HLTJpsi2Mu_L3_size);
   fChain->SetBranchAddress("HLTJpsi2Mu_L3_4mom", &HLTJpsi2Mu_L3_4mom, &b_HLTJpsi2Mu_L3_4mom);
   fChain->SetBranchAddress("HLTJpsi2Mu_L3_id", HLTJpsi2Mu_L3_id, &b_HLTJpsi2Mu_L3_id);
   fChain->SetBranchAddress("HLTUpsilon2Mu_L3_size", &HLTUpsilon2Mu_L3_size, &b_HLTUpsilon2Mu_L3_size);
   fChain->SetBranchAddress("HLTUpsilon2Mu_L3_4mom", &HLTUpsilon2Mu_L3_4mom, &b_HLTUpsilon2Mu_L3_4mom);
   fChain->SetBranchAddress("HLTUpsilon2Mu_L3_id", HLTUpsilon2Mu_L3_id, &b_HLTUpsilon2Mu_L3_id);
   Notify();
}

Bool_t JPsiTreeBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JPsiTreeBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JPsiTreeBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JPsiTreeBase_cxx
