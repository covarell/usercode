/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooRapiditySig.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "RooComplex.h"
#include "RooMath.h"
#include "TMath.h" 

ClassImp(RooRapiditySig) 

RooRapiditySig::RooRapiditySig(const char *name, const char *title, 
									   RooAbsReal& _Y,
			                                                   RooAbsReal& _m,
			                                                   RooAbsReal& _sqrtS) :
RooAbsPdf(name,title), 
Y("Y","Y",this,_Y),
m("m","m",this,_m),
sqrtS("sqrtS","sqrtS",this,_sqrtS)
{ 
} 


RooRapiditySig::RooRapiditySig(const RooRapiditySig& other, const char* name) :  
RooAbsPdf(other,name), 
Y("Y",this,other.Y),
m("m",this,other.m),
sqrtS("sqrtS",this,other.sqrtS)
{ 
} 



Double_t RooRapiditySig::evaluate() const 
{ 
        Double_t s0 = sqrtS*sqrtS;
	Double_t s = m*m;
	Double_t Q = m;
	Double_t xa = exp(Y)*sqrt(s/s0);
	Double_t xb = exp(-Y)*sqrt(s/s0);

	Double_t weightg = 1.0;

	//gluon params
	Double_t g0par0 = 0.2282; Double_t g0par1 = -0.0002252; Double_t g0par2 = 1.383e-07;
	Double_t g1par0 = 0.01968; Double_t g1par1 = -0.0002993; Double_t g1par2 = 1.986e-07;
	Double_t g2par0 = 3.624; Double_t g2par1 = -0.003164; Double_t g2par2 = 1.941e-06;
	Double_t g3par0 = -0.578; Double_t g3par1 = -0.0003; Double_t g3par2 = 1.828e-07;
	Double_t g4par0 = -7.515; Double_t g4par1 = -0.001355; Double_t g4par2 = 8.199e-07;

	Double_t gluon0 = g0par0 + g0par1*Q + g0par2*Q*Q;
	Double_t gluon1 = g1par0 + g1par1*Q + g1par2*Q*Q;
	Double_t gluon2 = g2par0 + g2par1*Q + g2par2*Q*Q;
	Double_t gluon3 = g3par0 + g3par1*Q + g3par2*Q*Q;
	Double_t gluon4 = g4par0 + g4par1*Q + g4par2*Q*Q;

	Double_t Funcga = (gluon0+gluon1*xa+gluon2*pow(xa,2))*pow((1-xa),4)*pow(xa,gluon3)*exp(1.0+gluon4*xa);
	Double_t Funcgb = (gluon0+gluon1*xb+gluon2*pow(xb,2))*pow((1-xb),4)*pow(xb,gluon3)*exp(1.0+gluon4*xb);
	Double_t FuncABg = Funcga*Funcgb/xa/xb;
	
	Double_t totSec = 2*m*((FuncABg)*weightg);

	if(( m <= 600. && TMath::Abs(Y) > 20*pow(m,-0.32)) || ( m > 600. && TMath::Abs(Y) > 21*pow(m,-0.34)))
	  {
	    //Find totSec when mZZ, Y=0
	    Double_t xa0 = sqrt(s/s0); //at Y=0 xa=xb
	    //if xa=xb then Funcga=Funcgb
	    Funcga = (gluon0+gluon1*xa0+gluon2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,gluon3)*exp(1.0+gluon4*xa0);
	    FuncABg = Funcga*Funcga/xa0/xa0;
	    Double_t totSec0 = 2*m*((FuncABg)*weightg);
	    totSec = 1.e-5*totSec0;
	  }

       	return totSec;	                                        
       
	
	
} 






