#ifndef MakeDataSet_h
#define MakeDataSet_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


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
  int get_Jpsi_pt_type(const double jpsi4mom);

private:

  int MIN_nhits_trk;
  float MAX_normchi2_trk;
  float MAX_normchi2_glb;
  int MIN_nhits_pixel;
  float MAX_d0_trk;
  float MAX_dz_trk;
  float MIN_vtxprob_jpsi;
  
};
#endif

