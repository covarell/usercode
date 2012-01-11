/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_PENTASPINTWO
#define ROO_PENTASPINTWO

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooPentaSpinTwo : public RooAbsPdf {
public:
  RooPentaSpinTwo(const char *name, const char *title,
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
		  RooAbsReal& _N);
  RooPentaSpinTwo(const RooPentaSpinTwo& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPentaSpinTwo(*this,newname); }
  inline virtual ~RooPentaSpinTwo() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy h1 ;
  RooRealProxy h2 ;
  RooRealProxy Phi ;
  RooRealProxy hs ;
  RooRealProxy Phi1 ;
  RooRealProxy fppVal ;
  RooRealProxy fmmVal ;
  RooRealProxy fpmVal ;
  RooRealProxy fp0Val ;
  RooRealProxy f0mVal ;
  RooRealProxy phippVal ;
  RooRealProxy phimmVal ;
  RooRealProxy phipmVal ;
  RooRealProxy phip0Val ;
  RooRealProxy phi0mVal ;
  RooRealProxy fz1Val ;
  RooRealProxy fz2Val ;
  RooRealProxy R1Val ;
  RooRealProxy R2Val ;
  RooRealProxy para2 ;
  RooRealProxy para4 ;
  RooRealProxy para6 ;
  RooRealProxy para8 ;
  RooRealProxy acca0 ;
  RooRealProxy acca1 ;
  RooRealProxy acca2 ;
  RooRealProxy acca4 ;
  RooRealProxy a2 ;
  RooRealProxy a4 ;
  RooRealProxy cutOff ;
  RooRealProxy g ;
  RooRealProxy b2 ;
  RooRealProxy b4 ;
  RooRealProxy N  ;
  Double_t evaluate() const ;

private:

  // ClassDef(RooPentaSpinTwo,1) // Your description goes here...
};
 
#endif
