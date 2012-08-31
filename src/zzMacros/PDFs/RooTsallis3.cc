
#include <iostream>
#include <math.h>

#include "TH1.h"
#include "RooFit.h"
#include "RooTsallis3.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

ClassImp(RooTsallis3) 

 RooTsallis3::RooTsallis3(const char *name, const char *title, 
				  RooAbsReal& _x,
				  RooAbsReal& _m,
				  RooAbsReal& _n,
                                  RooAbsReal& _n2,
                                  RooAbsReal& _bb,
			          RooAbsReal& _bb2,
			          RooAbsReal& _T,
                                  RooAbsReal& _fexp):
   
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   m("m","m",this,_m),
   n("n","n",this,_n),
   n2("n2","n2",this,_n2),
   bb("bb","bb",this,_bb),
   bb2("bb2","bb2",this,_bb2),
   T("T","T",this,_T),
   fexp("fexp","fexp",this,_fexp)
 { 
 } 


 RooTsallis3::RooTsallis3(const RooTsallis3& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   m("m",this,other.m),
   n("n",this,other.n),
   n2("n2",this,other.n2),
   bb("bb",this,other.bb),
   bb2("bb2",this,other.bb2),
   T("T",this,other.T),
   fexp("fexp",this,other.fexp)
 {
 } 



 double RooTsallis3::evaluate() const 
 { 
   // cout<<"In rooTsallis3::evaluate()"<<endl;
   return pow(x,n2)*exp(-bb*x)*pow(1 + (sqrt(x*x + m*m) - m)/(n*T),-n) + fexp*exp(-bb2*x);
 } 


// LET ROOFIT COMPUTE IT

/* int RooTsallis3::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}

double RooTsallis3::analyticalIntegral(int code, const char* rangeName) const
{
  switch(code)
    {
    case 1:
      {
	// mathematica dixit
	float term1 = x*pow(1 + (sqrt(x*x + m*m) - m)/(n*T),-n);
        float term2 = m*n*(sqrt(x*x + m*m) + 2*T) - n*n*T*(sqrt(x*x + m*m) + T) -n*m*m - n*x*x + x*x;
	return -term1*term2/((n-2)*(n-1));
      }
    }

assert(0) ;
return 0 ;
} */





