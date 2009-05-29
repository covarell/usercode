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
  
  // histos
  TH1F *QQMass2Glob_passmu3;
  TH1F *QQMass1Glob1Trk_passmu3;
  TH1F *QQMass1Glob1Cal_passmu3;
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
  
};
#endif

