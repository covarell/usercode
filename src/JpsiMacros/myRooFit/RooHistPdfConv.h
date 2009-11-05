/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooHistPdfConv.h,v 1.21 2007/05/11 09:13:07 verkerke Exp $
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
#ifndef ROO_MYGAUSS_MODEL
#define ROO_MYGAUSS_MODEL

#include "RooResolutionModel.h"
#include "RooRealProxy.h"
#include "RooDataHist.h"
#include "RooComplex.h"
#include "RooMath.h"

class RooHistPdfConv : public RooResolutionModel {
public:

  // enum RooGaussBasis { histBasis=1 };

  // Constructors, assignment etc
  inline RooHistPdfConv() { }
  RooHistPdfConv(const char *name, const char *title, RooRealVar& x, 
		RooAbsReal& mean, RooAbsReal& sigma, RooDataHist& datahist) ; 
  RooHistPdfConv(const char *name, const char *title, RooRealVar& x, 
		RooAbsReal& mean, RooAbsReal& sigma, RooAbsReal& msSF, RooDataHist& datahist ) ; 
  RooHistPdfConv(const char *name, const char *title, RooRealVar& x, 
		RooAbsReal& mean, RooAbsReal& sigma, RooAbsReal& meanSF, RooAbsReal& sigmaSF, RooDataHist& datahist) ; 
  RooHistPdfConv(const RooHistPdfConv& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooHistPdfConv(*this,newname) ; }
  virtual ~RooHistPdfConv();
  
  virtual Int_t basisCode(const char* name) const ;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

  void advertiseFlatScaleFactorIntegral(Bool_t flag) { _flatSFInt = flag ; }

  void advertiseAymptoticIntegral(Bool_t flag) { _asympInt = flag ; }  // added FMV,07/24/03

protected:

  virtual Double_t evaluate() const ;
  RooComplex evalCerfApprox(Double_t swt, Double_t u, Double_t c) const ;

  // Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
  inline RooComplex evalCerf(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFunc(z)*exp(-u*u) : evalCerfApprox(swt,u,c) ;
  }
    
  // Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  inline Double_t evalCerfRe(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncRe(z)*exp(-u*u) : evalCerfApprox(swt,u,c).re() ;
  }
  
  // Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  inline Double_t evalCerfIm(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncIm(z)*exp(-u*u) : evalCerfApprox(swt,u,c).im() ;
  }

  // Calculate Re(exp(-u^2) cwerf(i(u+c)))
  // added FMV, 08/17/03
  inline Double_t evalCerfRe(Double_t u, Double_t c) const {
    return exp(u*2*c+c*c) * RooMath::erfc(u+c);
  }

  // Calculate common normalization factors 
  // added FMV,07/24/03
  RooComplex evalCerfInt(Double_t sign, Double_t wt, Double_t tau, Double_t umin, Double_t umax, Double_t c) const ;
  Double_t evalCerfInt(Double_t sign, Double_t tau, Double_t umin, Double_t umax, Double_t c) const ;

  Bool_t _flatSFInt ;

  Bool_t _asympInt ;  // added FMV,07/24/03
  
  RooRealProxy mean ;
  RooRealProxy sigma ;
  RooRealProxy msf ;
  RooRealProxy ssf ;

  RooDataHist* _histpdf;
  const char* _variableName;

  ClassDef(RooHistPdfConv,1) // Gaussian Resolution Model
};

#endif









