
#include <iostream>
#include <math.h>

#include "RooPentaSpinTwo.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

// ClassImp(RooPentaSpinTwo) 

RooPentaSpinTwo::RooPentaSpinTwo(const char *name, const char *title, 
				  RooAbsReal& _h1,
				  RooAbsReal& _h2,
				  RooAbsReal& _Phi,
                                  RooAbsReal& _hs,
				  RooAbsReal& _Phi1,
				  RooAbsReal& _fppVal,
				  RooAbsReal& _fmmVal,
				  RooAbsReal& _fpmVal,
				  RooAbsReal& _fp0Val,
				  RooAbsReal& _f0mVal,
				  RooAbsReal& _phippVal,
				  RooAbsReal& _phimmVal,
				  RooAbsReal& _phipmVal,
				  RooAbsReal& _phip0Val,
				  RooAbsReal& _phi0mVal,
				  RooAbsReal& _fz1Val,
				  RooAbsReal& _fz2Val,
				  RooAbsReal& _R1Val,
				  RooAbsReal& _R2Val,
				  RooAbsReal& _para2,
				  RooAbsReal& _para4,
				  RooAbsReal& _para6,
				  RooAbsReal& _para8,
				  RooAbsReal& _acca0,
				  RooAbsReal& _acca1,
				  RooAbsReal& _acca2,
				  RooAbsReal& _acca4,
				  RooAbsReal& _a2,
				  RooAbsReal& _a4,
				  RooAbsReal& _cutOff,
				  RooAbsReal& _g,
                                  RooAbsReal& _b2,
				  RooAbsReal& _b4,
				  RooAbsReal& _N):
   
   RooAbsPdf(name,title), 
   h1("h1","h1",this,_h1),
   h2("h2","h2",this,_h2),
   Phi("Phi","Phi",this,_Phi),
   hs("hs","hs",this,_hs),
   Phi1("Phi1","Phi1",this,_Phi1),
   fppVal("fppVal","fppVal",this,_fppVal),
   fmmVal("fmmVal","fmmVal",this,_fmmVal),
   fpmVal("fpmVal","fpmVal",this,_fpmVal),
   fp0Val("fp0Val","fp0Val",this,_fp0Val),
   f0mVal("f0mVal","f0mVal",this,_f0mVal),
   phippVal("phippVal","phippVal",this,_phippVal),
   phimmVal("phimmVal","phimmVal",this,_phimmVal),
   phipmVal("phipmVal","phipmVal",this,_phipmVal),
   phip0Val("phip0Val","phip0Val",this,_phip0Val),
   phi0mVal("phi0mVal","phi0mVal",this,_phi0mVal),
   fz1Val("fz1Val","fz1Val",this,_fz1Val),
   fz2Val("fz2Val","fz2Val",this,_fz2Val),
   R1Val("R1Val","R1Val",this,_R1Val),
   R2Val("R2Val","R2Val",this,_R2Val),
   para2("para2","para2",this,_para2),
   para4("para4","para4",this,_para4),
   para6("para6","para6",this,_para6),
   para8("para8","para8",this,_para8),
   acca0("acca0","acca0",this,_acca0),
   acca1("acca1","acca1",this,_acca1),
   acca2("acca2","acca2",this,_acca2),
   acca4("acca4","acca4",this,_acca4),
   a2("a2","a2",this,_a2),
   a4("a4","a4",this,_a4),
   cutOff("cutOff","cutOff",this,_cutOff),
   g("g","g",this,_g),
   b2("b2","b2",this,_b2),
   b4("b4","b4",this,_b4),
   N("N","N",this,_N)
 { 
 } 


 RooPentaSpinTwo::RooPentaSpinTwo(const RooPentaSpinTwo& other, const char* name) :  
   RooAbsPdf(other,name), 
   h1("h1",this,other.h1),
   h2("h2",this,other.h2),
   Phi("Phi",this,other.Phi),
   hs("hs",this,other.hs),
   Phi1("Phi1",this,other.Phi1),
   fppVal("fppVal",this,other.fppVal),
   fmmVal("fmmVal",this,other.fmmVal),
   fpmVal("fpmVal",this,other.fpmVal),
   fp0Val("fp0Val",this,other.fp0Val),
   f0mVal("f0mVal",this,other.f0mVal),
   phippVal("phippVal",this,other.phippVal),
   phimmVal("phimmVal",this,other.phimmVal),
   phipmVal("phipmVal",this,other.phipmVal),
   phip0Val("phip0Val",this,other.phip0Val),
   phi0mVal("phi0mVal",this,other.phi0mVal),
   fz1Val("fz1Val",this,other.fz1Val),
   fz2Val("fz2Val",this,other.fz2Val),
   R1Val("R1Val",this,other.R1Val),
   R2Val("R2Val",this,other.R2Val),
   para2("para2",this,other.para2),
   para4("para4",this,other.para4),
   para6("para6",this,other.para6),
   para8("para8",this,other.para8),
   acca0("acca0",this,other.acca0),
   acca1("acca1",this,other.acca1),
   acca2("acca2",this,other.acca2),
   acca4("acca4",this,other.acca4),
   a2("a2",this,other.a2),
   a4("a1",this,other.a4),
   cutOff("cutOff",this,other.cutOff),
   g("g",this,other.g),
   b2("b2",this,other.b2),
   b4("b4",this,other.b4),
   N("N",this,other.N)
 {
 } 



 double RooPentaSpinTwo::evaluate() const 
 { 
   // cout<<"In rooPentaSpinTwo::evaluate()"<<endl;
   double shs = sqrt(1-hs*hs);
   double sh1 = sqrt(1-h1*h1);
   double sh2 = sqrt(1-h2*h2);
   
   if ((1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal) < 0) return 1e-9;
   
   double term1Coeff = (2.-2.*fz1Val+fz2Val-6.*(2.-4.*fz1Val-fz2Val)*pow(hs,2)+3.*(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term1A = 4.*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal)*pow(sh1,2)*pow(sh2,2);
   double term1B = (fppVal+fmmVal)*((1.+h1*h1)*(1.+h2*h2)+4.*R1Val*R2Val*h1*h2);
   double term1C = -2.*(fppVal-fmmVal)*(R1Val*h1*(1.+h2*h2)+R2Val*h2*(1.+h1*h1));
   double term1D = 4.*sqrt(fppVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val-h1)*(R2Val-h2)*sh1*sh2*cos(Phi+phippVal);
   double term1E = 4.*sqrt(fmmVal*(1-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(R1Val+h1)*(R2Val+h2)*sh1*sh2*cos(Phi-phimmVal);
   double term1F = 2.*sqrt(fppVal*fmmVal)*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+phippVal-phimmVal);
   double term1 = term1Coeff*(term1A+term1B+term1C+term1D+term1E+term1F);
   // cout<<"CP1  "<<flush;   
   double term2Coeff = 8.*(fz1Val+fz2Val+3.*(2.-3.*fz1Val-2.*fz2Val)*pow(hs,2)-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4));
   double term2A = (fp0Val+f0mVal)*(1.-h1*h1*h2*h2) - (fp0Val-f0mVal)*(R1Val*h1*pow(sh2,2)+R2Val*h2*pow(sh1,2));
   double term2B = 2.*sqrt(fp0Val*f0mVal)*sh1*sh2*(R1Val*R2Val-h1*h2)*cos(Phi+phip0Val-phi0mVal);
   double term2 = term2Coeff*(term2A+term2B);
   
   double term3Coeff = -8.*(fz1Val-fz2Val+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2)*sh1*sh2*cos(Phi + 2.*Phi1);
   double term3A = (fp0Val+f0mVal)*(R1Val*R2Val+h1*h2)-(fp0Val-f0mVal)*(R1Val*h2+R2Val*h1)+2.*sqrt(fp0Val*f0mVal)*sh1*sh2*cos(Phi+phip0Val-phi0mVal);
   double term3 = term3Coeff*term3A;
   
   double term4Coeff = 6.-2.*fz1Val-5.*fz2Val-6.*(2.-2.*fz1Val-3.*fz2Val)*pow(hs,2)+(6.-10.*fz1Val-5.*fz2Val)*pow(hs,4);
   double term4A = fpmVal*((1.+h1*h1)*(1.+h2*h2)-4.*R1Val*R2Val*h1*h2);
   double term4 = term4Coeff*term4A;
    // cout<<"CP2  "<<flush;   
   double term5 = pow(shs,4)*(6.-10.*fz1Val-5.*fz2Val)*fpmVal*pow(sh1,2)*pow(sh2,2)*cos(2.*Phi+4.*Phi1);
   
   double term6Coeff = (-1.)*sqrt(6.)*(2.-2.*fz1Val-3.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*pow(shs,2);
   double term6A = 2.*sqrt(fpmVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*sh1*sh2*((R1Val-h1)*(R2Val+h2)*cos(Phi+2.*Phi1-phipmVal)+(R1Val+h1)*(R2Val-h2)*cos(Phi+2.*Phi1+phipmVal));
   double term6B = sqrt(fpmVal*fmmVal)*(sh1*sh1*(1.+2.*R2Val*h2+h2*h2)*cos(2.*Phi1-phipmVal+phimmVal)+sh2*sh2*(1.+2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1+phipmVal-phimmVal));
   double term6C = sqrt(fppVal*fpmVal)*(sh1*sh1*(1.-2.*R2Val*h2+h2*h2)*cos(2.*Phi1+phipmVal-phippVal)+sh2*sh2*(1.-2.*R1Val*h1+h1*h1)*cos(2.*Phi+2.*Phi1-phipmVal+phippVal));
   double term6 = term6Coeff*(term6A+term6B+term6C);
   
   /// new mixing terms
   double term7Coeff = -4.*sqrt(3.)*(2.-4.*fz1Val-fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2))*hs*shs;
   double term7A = sqrt(fmmVal*f0mVal)*(sh1*(R1Val+h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1+phimmVal-phi0mVal)-sh2*(R2Val+h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phimmVal+phi0mVal));
   double term7B = sqrt(fmmVal*fp0Val)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi-Phi1-phimmVal+phip0Val)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+Phi1-phimmVal+phip0Val));
   double term7C = -sqrt(fppVal*f0mVal)*(sh1*sh1*sh2*(R2Val-h2)*cos(Phi-Phi1+phippVal-phi0mVal)-sh2*sh2*sh1*(R1Val-h1)*cos(2.*Phi+Phi1+phippVal-phi0mVal));
   double term7D = -sqrt(fppVal*fp0Val)*(sh1*(R1Val-h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phip0Val-phippVal)-sh2*(R2Val-h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phippVal-phip0Val));
   double term7E = -2.*sqrt(f0mVal*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val+h1)*cos(Phi1+phi0mVal)-sh2*sh1*sh1*(R2Val+h2)*cos(Phi+Phi1-phi0mVal));
   double term7F = 2.*sqrt(fp0Val*(1.-fppVal-fmmVal-2.*fpmVal-2.*fp0Val-2.*f0mVal))*(sh1*sh2*sh2*(R1Val-h1)*cos(Phi1-phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+Phi1+phip0Val));
   double term7 = term7Coeff*(term7A+term7B+term7C+term7D+term7E+term7F);
    // cout<<"CP3  "<<flush;   
   double term8Coeff = 2.*sqrt(2.)*hs*shs*(6.-6.*fz1Val-9.*fz2Val-(6.-10.*fz1Val-5.*fz2Val)*pow(hs,2));
   double term8A = sqrt(fpmVal*f0mVal)*(sh1*(R1Val-h1)*(1.+2.*R2Val*h2+h2*h2)*cos(Phi1-phipmVal+phi0mVal)-sh2*(R2Val-h2)*(1.+2.*R1Val*h1+h1*h1)*cos(Phi+Phi1+phipmVal-phi0mVal));
   double term8B = sqrt(fpmVal*fp0Val)*(sh2*(R2Val+h2)*(1.-2.*R1Val*h1+h1*h1)*cos(Phi+Phi1-phipmVal+phip0Val)-sh1*(R1Val+h1)*(1.-2.*R2Val*h2+h2*h2)*cos(Phi1+phipmVal-phip0Val));
   double term8 = term8Coeff*(term8A+term8B);
   
   double term9Coeff = -2.*sqrt(2.)*hs*pow(shs,3)*(6.-10.*fz1Val-5.*fz2Val);
   double term9A = sqrt(fpmVal*f0mVal)*(sh1*sh1*sh2*(R2Val+h2)*cos(Phi+3.*Phi1-phipmVal+phi0mVal)-sh2*sh2*sh1*(R1Val+h1)*cos(2.*Phi+3.*Phi1+phipmVal-phi0mVal));
   double term9B = sqrt(fpmVal*fp0Val)*(sh1*sh2*sh2*(R1Val-h1)*cos(2.*Phi+3.*Phi1-phipmVal+phip0Val)-sh2*sh1*sh1*(R2Val-h2)*cos(Phi+3.*Phi1+phipmVal-phip0Val));
   double term9 = term9Coeff*(term9A+term9B);
    // cout<<"CP4  "<<flush;   
   // signs of the interference terms are flipped!!!!
   if (true){
     term7 *= (-1.);
     term8 *= (-1.);
     term9 *= (-1.);
   }
   
   double sum = term1+term2+term3+term4+term5+term6+term7+term8+term9;

   double accp = (1.0+para2*pow(hs,2)+para4*pow(hs,4)+para6*pow(hs,6)+para8*pow(hs,8))
                  *(1+b2*pow(h1,2)+b4*pow(h1,4))
                  *(acca0+acca1*cos(Phi1)+acca2*cos(2.0*Phi1)+acca4*cos(4.0*Phi1))
                  *(1.0+a2*h2+a4*pow(h2,2))/(1.0+exp((h2-cutOff)*g));
                  // TEST ANOTHER SHAPE (cutOff = a1, g = a3)
                  // *(1.0+cutOff*h2+a2*pow(h2,2)+g*pow(h2,3)+a4*pow(h2,4));

   if(sum*accp<=0){sum=accp=.0001;}
   // cout<<"after all magics: h2="<<h2<<"  accp="<<accp<<endl;
   return sum*accp ;
 } 

int RooPentaSpinTwo::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars,analVars,h1,Phi1,Phi)) return 1 ;
  return 0 ;
}

double RooPentaSpinTwo::analyticalIntegral(int code, const char* rangeName) const
{
  switch(code)
    {
    case 1:
      {
	double phi00Val=0;
	double f00;
	f00=1-fmmVal-fppVal-2*f0mVal-2*fp0Val-2*fpmVal;
	double term1 =  (-2*acca0*3.1415*(4*(35 + 7*b2 + 3*b4)*f00*3.1415*(-1 + pow(h2,2)) - 2*(35 + 14*b2 + 9*b4)*3.1415*(fmmVal + fppVal + 2*(fmmVal - fppVal)*R2Val*h2 + (fmmVal + fppVal)*pow(h2,2)))*(2 - 2*fz1Val + fz2Val + 6*(-2 + 4*fz1Val + fz2Val)*pow(hs,2) - 3*(-6 + 10*fz1Val + 5*fz2Val)*pow(hs,4)))/105;
	  double term2 = (8*acca0*3.1415*3.1415*(7*(5*b2 + 3*(5 + b4))*(f0mVal + fp0Val) + 2*(35 + 7*b2 + 3*b4)*(f0mVal - fp0Val)*R2Val*h2 - (35 + 21*b2 + 15*b4)*(f0mVal + fp0Val)*pow(h2,2))*(fz1Val + fz2Val - 3*(-2 + 3*fz1Val + 2*fz2Val)*pow(hs,2) + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,4)))/105;
	  double term3 = (-8*acca2*(35 + 7*b2 + 3*b4)*sqrt(f0mVal)*sqrt(fp0Val)*3.1415*3.1415*(1 - pow(h2,2))*(-1 + pow(hs,2))*(-fz1Val + fz2Val + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,2))*cos(phi0mVal - phip0Val))/105;
	  double term4 = (4*acca0*(35 + 14*b2 + 9*b4)*fpmVal*3.1415*3.1415*(1 + pow(h2,2))*(6 - 2*fz1Val - 5*fz2Val + 6*(-2 + 2*fz1Val + 3*fz2Val)*pow(hs,2) + (6 - 10*fz1Val - 5*fz2Val)*pow(hs,4)))/105;
	  double term5 = (sqrt(2/3)*acca2*(35 + 7*b2 + 3*b4)*sqrt(fpmVal)*3.1415*3.1415*(-1 + pow(hs,2))*(2 - 6*pow(hs,2) + fz2Val*(-3 + 5*pow(hs,2)) + 2*fz1Val*(-1 + 5*pow(hs,2)))*(sqrt(fmmVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phimmVal - phipmVal) + sqrt(fppVal)*(1 - 2*R2Val*h2 + pow(h2,2))*cos(phipmVal - phippVal)))/35;
	  double term6 = -(acca1*(8 + 2*b2 + b4)*3.1415*3.1415*3.1415*R1Val*hs*sqrt(3 - 3*pow(hs,2))*(-2 + 4*fz1Val + fz2Val + (6 - 10*fz1Val - 5*fz2Val)*pow(hs,2))*(2*sqrt(f00*f0mVal)*(-1 + pow(h2,2))*cos(phi00Val - phi0mVal) + sqrt(f0mVal)*sqrt(fmmVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phi0mVal - phimmVal) + 2*sqrt(f00*fp0Val)*cos(phi00Val - phip0Val) - 2*sqrt(f00*fp0Val)*pow(h2,2)*cos(phi00Val - phip0Val) - sqrt(fp0Val)*sqrt(fppVal)*cos(phip0Val - phippVal) + 2*sqrt(fp0Val)*sqrt(fppVal)*R2Val*h2*cos(phip0Val - phippVal) - sqrt(fp0Val)*sqrt(fppVal)*pow(h2,2)*cos(phip0Val - phippVal)))/16;
	  double term7 = -(acca1*(8 + 2*b2 + b4)*sqrt(fpmVal)*3.1415*3.1415*3.1415*R1Val*hs*sqrt(1 - pow(hs,2))*(6 - 6*fz1Val - 9*fz2Val + (-6 + 10*fz1Val + 5*fz2Val)*pow(hs,2))*(sqrt(f0mVal)*(1 + 2*R2Val*h2 + pow(h2,2))*cos(phi0mVal - phipmVal) + sqrt(fp0Val)*(-1 + 2*R2Val*h2 - pow(h2,2))*cos(phip0Val - phipmVal)))/(16*sqrt(2));

	  return N*(term1+term2+term3+term4+term5+term6+term7)*(1+para2*hs*hs+para4*hs*hs*hs*hs+para6*pow(hs,6)+para8*pow(hs,8))
	    *(1+a2*h2+a4*h2*h2)/(1 + exp((h2-cutOff)*g));
	    // TEST ANOTHER SHAPE (cutOff = a1, g = a3)
	    // *(1+cutOff*h2+a2*h2*h2+g*h2*h2*h2+a4*h2*h2*h2*h2);
      }
    }

assert(0) ;
return 0 ;
}





