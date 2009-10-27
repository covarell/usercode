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

private:

  int MIN_nhits_trk;
  float MAX_normchi2_trk;
  float MAX_normchi2_glb;
  int MIN_nhits_pixel;
  float MAX_d0_trk;
  float MAX_dz_trk;
  
};
#endif

