#include <vector>
#include "TLorentzVector.h"
#include "RooWorkspace.h"

// void initAllFunctions(RooWorkspace *ws, int LHCsqrts);

template <typename U>
void checkZorder(U& z1mass, U& z2mass,
                 U& costhetastar, U& costheta1,
                 U& costheta2, U& phi,
                 U& phistar1, U& pt, U& ipsilon);

vector<double> my8DTemplate(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1);

pair<double,double> likelihoodDiscriminant (RooWorkspace *ws,
					    double mZZ, double m1, double m2, 
					    double costhetastar, double costheta1, double costheta2, double phi, double phi1, 
					    int LHCsqrts, 
					    bool withPt = false, double pt = 0.0, 
					    bool withY = false, double absy = 0.0, 
					    double scaleFactor=5.0);

void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4Lep11, TLorentzVector thep4Lep12, TLorentzVector thep4Z2, TLorentzVector thep4Lep21, TLorentzVector thep4Lep22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2);

void setTemplate(char* file);
