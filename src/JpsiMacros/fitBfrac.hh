#ifndef fitBfrac_h
#define fitBfrac_h

#include <TProfile2D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <map>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "JPsiTreeBase.h"


class fitBfrac : public JPsiTreeBase{
public:
  
  fitBfrac(TTree *tree=0);
  virtual ~fitBfrac();
  void Loop();

private:

  void bookHistos();
  void drawPlots();
  void saveHistos();
  
  // histos
  // TH1F *QQMass2Glob_passmu3;
  
  
};
#endif

