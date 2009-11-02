#ifndef MakeDataSet_h
#define MakeDataSet_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"
#include <vector>

class MakeDataSet : public JPsiTreeBase{
public:
  
  MakeDataSet(TTree *tree=0);
  virtual ~MakeDataSet();
  void Loop();
  int theBestQQ();
  double PhiInRange(const double& phi) const;
  double deltaR(const TLorentzVector* t, const TLorentzVector* u) const;
  float findEff(TH2F* effhist, float pt, float eta, bool approx) const;
  float findEffErr(TH2F* effhist, float pt, float eta, bool approx) const;
  int get_Jpsi_var_type(const double jpsi4mom, vector<double> min, vector<double> max);

private:

  int MIN_nhits_trk;
  float MAX_normchi2_trk;
  float MAX_normchi2_glb;
  int MIN_nhits_pixel;
  float MAX_d0_trk;
  float MAX_dz_trk;
  float MIN_vtxprob_jpsi;

  TH1F* gammaFactor_GGnonprompt;
  TH1F* gammaFactor_GTnonprompt;
  TH1F* hMcPR_GGMass;               
  TH1F* hMcPR_GTMass;               
  TH1F* hMcNP_GGMass;               
  TH1F* hMcNP_GTMass;               
  TH1F* hMcBK_GGMass;               
  TH1F* hMcBK_GTMass;               
  TH1F* hMcPR_GGLife;               
  TH1F* hMcPR_GTLife;               
  TH1F* hMcNP_GGLife;               
  TH1F* hMcNP_GTLife;               
  TH1F* hMcBK_GGLife;               
  TH1F* hMcBK_GTLife; 

  TH2F *heffTrk;
  TH2F *heffMuGlb;
  TH2F *heffMuTrk;
  TH2F *heffMuHLT;

  bool onlyTheBest;
  bool efficiencyStore;
              
};
#endif

