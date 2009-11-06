/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooHistPdfConv.cxx,v 1.1 2009/11/05 16:38:57 covarell Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooHistPdfConv implements a RooResolutionModel that models a Gaussian
// distribution. Object of class RooHistPdfConv can be used
// for analytical convolutions with classes inheriting from RooAbsAnaConvPdf
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include "RooHistPdfConv.h"
#include "RooMath.h"
#include "RooRealConstant.h"
#include "RooRandom.h"
#include "RooRealVar.h"

ClassImp(RooHistPdfConv);


//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooDataHist& datahist) :
  RooAbsPdf(name,title), 
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  xIn("xIn","xIn",this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,(RooAbsReal&)RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,(RooAbsReal&)RooRealConstant::value(1))
{  
  _histpdf = new RooDataHist(datahist);
  _variableName = xIn.GetName();
}



//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _msSF, RooDataHist& datahist) : 
  RooAbsPdf(name,title), 
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  xIn("xIn","xIn",this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_msSF),
  ssf("ssf","Sigma Scale Factor",this,_msSF)
{
  _histpdf = new RooDataHist(datahist);
  _variableName = xIn.GetName();
}



//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _meanSF, RooAbsReal& _sigmaSF,
                             RooDataHist& datahist ) : 
  RooAbsPdf(name,title), 
  _flatSFInt(kFALSE),
  _asympInt(kFALSE),
  xIn("xIn","xIn",this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{   
  _histpdf = new RooDataHist(datahist); 
  _variableName = xIn.GetName();
}   


//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const RooHistPdfConv& other, const char* name) : 
  RooAbsPdf(other,name),
  _flatSFInt(other._flatSFInt),
  _asympInt(other._asympInt),
  xIn("xIn",this,other.xIn),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  msf("msf",this,other.msf),
  ssf("ssf",this,other.ssf)
{
  _histpdf = other._histpdf;
}



//_____________________________________________________________________________
RooHistPdfConv::~RooHistPdfConv()
{
  // Destructor
}


//_____________________________________________________________________________
Double_t RooHistPdfConv::evaluate() const 
{  
  // cout << "RooHistPdfConv::evaluate(" << GetName() << ")" << endl ;
  
  static Double_t root2(sqrt(2.)) ; 
  const RooArgSet* aRow;
  RooRealVar* xprime;
 
  Double_t result(0) ;
  // *** Convolution with hist PDF ***
  for (Int_t i=0; i<_histpdf->numEntries(); i++) {
    
    aRow = _histpdf->get(i);
    Double_t halfBinSize = _histpdf->binVolume(*aRow)/2.0;
    Double_t weight = _histpdf->weight(*aRow,0,false)/_histpdf->sum(false);
    xprime = (RooRealVar*)aRow->find(_variableName);

    Double_t c = (xprime->getVal() - halfBinSize - xIn + (mean*msf)) / (root2*sigma*ssf);
    Double_t d = (xprime->getVal() + halfBinSize - xIn + (mean*msf)) / (root2*sigma*ssf);
    result += 0.5*weight*(RooMath::erfc(c)-RooMath::erfc(d));
  }

  return result ;
}

//_____________________________________________________________________________
Int_t RooHistPdfConv::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,xIn)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t RooHistPdfConv::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  /* static Double_t root2 = sqrt(2.) ;
  //static Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
 
  // *** 3rd form: Convolution with exp(-t/tau), used for expBasis and cosBasis(omega=0) ***
  Double_t c = (sigma*ssf)/(root2*tau) ; 
  Double_t xpmin = (x.min(rangeName)-(mean*msf))/tau ;
  Double_t xpmax = (x.max(rangeName)-(mean*msf))/tau ;
  Double_t umin = xpmin/(2*c) ;
  Double_t umax = xpmax/(2*c) ;

  if (basisType==expBasis || (basisType==cosBasis && omega==0.)) {
    if (verboseEval()>0) cout << "RooHistPdfConv::analyticalIntegral(" << GetName() << ") 3d form tau=" << tau << endl ;

    Double_t result(0) ;
    if (_asympInt) {   // modified FMV, 07/24/03
      if (basisSign!=Minus) result += 2 * tau ;
      if (basisSign!=Plus)  result += 2 * tau ;      
    } else {
      if (basisSign!=Minus) result += -1 * tau * ( RooMath::erf(-umax) - RooMath::erf(-umin) + 
						   exp(c*c) * ( exp(-xpmax)*RooMath::erfc(-umax+c)
								- exp(-xpmin)*RooMath::erfc(-umin+c) )) ;
      if (basisSign!=Plus)  result +=      tau * ( RooMath::erf(umax) - RooMath::erf(umin) + 
						   exp(c*c) * ( exp(xpmax)*RooMath::erfc(umax+c)
								- exp(xpmin)*RooMath::erfc(umin+c) )) ;     
      // equivalent form, added FMV, 07/24/03
      //if (basisSign!=Minus) result += evalCerfInt(+1,tau,-umin,-umax,c).re();   
      //if (basisSign!=Plus) result += evalCerfInt(-1,tau,umin,umax,c).re();
    }
    //cout << "Integral 3rd form " << " result= " << result*ssfInt << endl;
    return result*ssfInt ; */
  
 
  assert(0) ;
  return 0 ;
}



//_____________________________________________________________________________
RooComplex RooHistPdfConv::evalCerfApprox(Double_t swt, Double_t u, Double_t c) const
{
  // use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z)
  // to explicitly cancel the divergent exp(y*y) behaviour of
  // CWERF for z = x + i y with large negative y

  static Double_t rootpi= sqrt(atan2(0.,-1.));
  RooComplex z(swt*c,u+c);  
  RooComplex zc(u+c,-swt*c);
  RooComplex zsq= z*z;
  RooComplex v= -zsq - u*u;

  return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ;
}



// added FMV, 07/24/03
//_____________________________________________________________________________
RooComplex RooHistPdfConv::evalCerfInt(Double_t sign, Double_t wt, Double_t tau, Double_t umin, Double_t umax, Double_t c) const
{
  RooComplex diff;
  if (_asympInt) {
    diff = RooComplex(2,0) ;
  } else {
    diff = RooComplex(sign,0.)*(evalCerf(wt,umin,c) - evalCerf(wt,umax,c) + RooMath::erf(umin) - RooMath::erf(umax));
  }
  return RooComplex(tau/(1.+wt*wt),0)*RooComplex(1,wt)*diff;
}
// added FMV, 08/17/03

//_____________________________________________________________________________
Double_t RooHistPdfConv::evalCerfInt(Double_t sign, Double_t tau, Double_t umin, Double_t umax, Double_t c) const
{
  Double_t diff;
  if (_asympInt) {
    diff = 2. ;
  } else {
    if ((umin<-8 && umax>8)||(umax<-8 && umin>8)) {
      // If integral is over >8 sigma, approximate with full integral
      diff = 2. ;
    } else {
      diff = sign*(evalCerfRe(umin,c) - evalCerfRe(umax,c) + RooMath::erf(umin) - RooMath::erf(umax));
    }
  }
  return tau*diff;
}



//_____________________________________________________________________________
Int_t RooHistPdfConv::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,xIn)) return 1 ;  
  return 0 ;
}



//_____________________________________________________________________________
void RooHistPdfConv::generateEvent(Int_t code)
{
  assert(code==1) ;
  Double_t xgen ;
  while(1) {
    xgen = RooRandom::randomGenerator()->Gaus((mean*msf),(sigma*ssf));
    if (xgen < xIn.max() && xgen > xIn.min()) {
      xIn = xgen ;
      return ;
    }
  }
}




