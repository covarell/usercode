#ifndef finalJPsiAnalysis_h
#define finalJPsiAnalysis_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


class finalJPsiAnalysis : public JPsiTreeBase{
public:
  
  finalJPsiAnalysis(TTree *tree=0);
  virtual ~finalJPsiAnalysis();
  void Loop();
  

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  int theBestQQ();
  double PhiInRange(const double& phi) const;
  double deltaR(const TLorentzVector* t, const TLorentzVector* u) const;

  // histos
  TH1F *QQMass2Glob_passmu3;
  TH1F *QQMass1Glob1Trk_passmu3;
  TH1F *QQMass1Glob1Cal_passmu3;
  TH1F *QQMass2Glob_best;
  TH1F *QQMass1Glob1Trk_best;
  TH1F *QQMass1Glob1Cal_best;
  TH1F *QQMass2GlobPT6_passmu3;
  TH1F *QQMass1Glob1TrkPT6_passmu3;
  TH1F *QQMass1Glob1CalPT6_passmu3;
  TH1F *WSMass2Glob_passmu3;
  TH1F *WSMass1Glob1Trk_passmu3;
  TH1F *WSMass1Glob1Cal_passmu3;
  TH1F *QQMass2Glob_passmu5;
  TH1F *QQMass1Glob1Trk_passmu5;
  TH1F *QQMass1Glob1Cal_passmu5;
  TH1F *QQMass2Glob_passmu9;
  TH1F *QQMass1Glob1Trk_passmu9;
  TH1F *QQMass1Glob1Cal_passmu9;
  TH1F *QQMass2GlobPT6_passmu5;
  TH1F *QQMass1Glob1TrkPT6_passmu5;
  TH1F *QQMass1Glob1CalPT6_passmu5;
  TH1F *QQMass2Glob_pass2mu3;
  TH1F *QQMass1Glob1Trk_pass2mu3;
  TH1F *QQMass1Glob1Cal_pass2mu3;
  TH1F *WSMass2Glob_pass2mu3;
  TH1F *WSMass1Glob1Trk_pass2mu3;
  TH1F *WSMass1Glob1Cal_pass2mu3;
  TH1F *QQMass2Glob_pass2isomu3;             
  TH1F *QQMass1Glob1Trk_pass2isomu3; 
  TH1F *QQMass1Glob1Cal_pass2isomu3;
  TH1F *QQMass2Glob_pass2mu3jpsi;    
  TH1F *QQMass1Glob1Trk_pass2mu3jpsi;
  TH1F *QQMass1Glob1Cal_pass2mu3jpsi;
  TH1F *QQMass2Glob_pass2mu3ups;     
  TH1F *QQMass1Glob1Trk_pass2mu3ups; 
  TH1F *QQMass1Glob1Cal_pass2mu3ups; 

  TH1F *hMcRecoGlobMuDeltaR;        
  TH1F *hMcRecoTrkMuDeltaR;  
  TH1F *hMcRecoCalMuDeltaR;  

  TH1F *hMcRightTrkMuPt;
  TH1F *hMcWrongTrkMuPt;
  TH1F *hMcRightTrkBit4;
  TH1F *hMcWrongTrkBit4;
  TH1F *hMcRightTrkBit5;
  TH1F *hMcWrongTrkBit5;
  TH1F *hMcRightTrkBit8;
  TH1F *hMcWrongTrkBit8;
  TH1F *hMcRightTrkBit9;
  TH1F *hMcWrongTrkBit9;
  TH1F *hMcRightTrkMuChi2;   
  TH1F *hMcWrongTrkMuChi2; 
  TH1F *hMcRightTrkMuNhits;
  TH1F *hMcWrongTrkMuNhits;

  TH1F *hMcRightCalMuChi2; 
  TH1F *hMcWrongCalMuChi2; 
  TH1F *hMcRightCalMuNhits;
  TH1F *hMcWrongCalMuNhits;
  TH1F *hMcRightCalMuCaloComp;         
  TH1F *hMcWrongCalMuCaloComp;      
  TH1F *hMcRightCalMuPt;            
  TH1F *hMcWrongCalMuPt;    
  TH1F *hMcRightCalGlobMuDeltaR;        
  TH1F *hMcWrongCalGlobMuDeltaR;     
  TH1F *hMcRightCalGlobMuMass;              
  TH1F *hMcWrongCalGlobMuMass;
  TH1F *hMcRightCalGlobMuVtxChi2; 
  TH1F *hMcWrongCalGlobMuVtxChi2; 
  TH1F *hMcRightCalGlobMuS; 
  TH1F *hMcWrongCalGlobMuS;
  TH1F *hMcRightCalGlobMucosAlpha; 
  TH1F *hMcWrongCalGlobMucosAlpha;
      
  TH1F *hMcRightAllMuIso;

};
#endif

