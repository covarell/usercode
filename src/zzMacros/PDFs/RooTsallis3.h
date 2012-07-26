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

#ifndef ROO_TSALLIS3
#define ROO_TSALLIS3

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooTsallis3 : public RooAbsPdf {
public:
  RooTsallis3(const char *name, const char *title,
	          RooAbsReal& _x,
        	  RooAbsReal& _m,
              	  RooAbsReal& _n,
	          RooAbsReal& _n2,
                  RooAbsReal& _bb,
                  RooAbsReal& _T);

  RooTsallis3(const RooTsallis3& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTsallis3(*this,newname); }
  inline virtual ~RooTsallis3() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x ;
  RooRealProxy m ;
  RooRealProxy n ;
  RooRealProxy n2 ;
  RooRealProxy bb ;
  RooRealProxy T ;
 
  Double_t evaluate() const ;

private:

 ClassDef(RooTsallis3,1) // Your description goes here...
};
 
#endif
