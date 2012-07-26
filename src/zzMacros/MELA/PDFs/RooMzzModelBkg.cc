#include <iostream>
#include <math.h>
#include <TMath.h>

#include "RooMzzModelBkg.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

 ClassImp(RooMzzModelBkg) 

 RooMzzModelBkg::RooMzzModelBkg(){}

 RooMzzModelBkg::RooMzzModelBkg(const char *name, const char *title,
				RooAbsReal& _mZZ, 
				RooAbsReal& _a0,
				RooAbsReal& _a1,
				RooAbsReal& _a2,
				RooAbsReal& _a3,
				RooAbsReal& _a4,
				RooAbsReal& _a5,
				RooAbsReal& _a6,
				RooAbsReal& _a7,
				RooAbsReal& _a8,
				RooAbsReal& _a9) :
   RooAbsPdf(name,title), 
   mZZ("mZZ","mZZ",this,_mZZ),
   a0("a0","a0",this,_a0),
   a1("a1","a1",this,_a1),
   a2("a2","a2",this,_a2),
   a3("a3","a3",this,_a3),
   a4("a4","a4",this,_a4),
   a5("a5","a5",this,_a5),
   a6("a6","a6",this,_a6),
   a7("a7","a7",this,_a7),
   a8("a8","a8",this,_a8),
   a9("a9","a9",this,_a9)
 { 

  } 

RooMzzModelBkg::RooMzzModelBkg(const RooMzzModelBkg& other, const char* name) :  
   RooAbsPdf(other,name), 
   mZZ("mZZ",this,other.mZZ),
   a0("a0",this,other.a0),
   a1("a1",this,other.a1),
   a2("a2",this,other.a2),
   a3("a3",this,other.a3),
   a4("a4",this,other.a4),
   a5("a5",this,other.a5),
   a6("a6",this,other.a6),
   a7("a7",this,other.a7),
   a8("a8",this,other.a8),
   a9("a9",this,other.a9)
 { 
   
 }

double RooMzzModelBkg::evaluate() const 
 { 
   double ZZ = (.5+.5*TMath::Erf((mZZ-a0)/a1))*(a3/(1+exp((mZZ-a0)/a2)))
     +(.5+.5*TMath::Erf((mZZ-a4)/a5))*(a7/(1+exp((mZZ-a4)/a6))+a9/(1+exp((mZZ-a4)/a8)));
     //(.5+.5*TMath::Erf((mZZ-a10)/a11))*(a13/(1+exp((mZZ-a10)/a12)) );
	
	return ZZ;
 } 

