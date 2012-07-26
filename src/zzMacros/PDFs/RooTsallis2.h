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

#ifndef ROO_TSALLIS2
#define ROO_TSALLIS2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooTsallis2 : public RooAbsPdf {
public:
  RooTsallis2(const char *name, const char *title,
	          RooAbsReal& _x,
        	  RooAbsReal& _m,
              	  RooAbsReal& _n,
	          RooAbsReal& _n2,
                  RooAbsReal& _T);

  RooTsallis2(const RooTsallis2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTsallis2(*this,newname); }
  inline virtual ~RooTsallis2() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x ;
  RooRealProxy m ;
  RooRealProxy n ;
  RooRealProxy n2 ;
  RooRealProxy T ;
 
  Double_t evaluate() const ;

private:

 ClassDef(RooTsallis2,1) // Your description goes here...
};
 
#endif
