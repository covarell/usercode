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

#ifndef ROO_TSALLIS
#define ROO_TSALLIS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooTsallis : public RooAbsPdf {
public:
  RooTsallis(const char *name, const char *title,
	          RooAbsReal& _x,
        	  RooAbsReal& _m,
              	  RooAbsReal& _n,
                  RooAbsReal& _T);

  RooTsallis(const RooTsallis& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTsallis(*this,newname); }
  inline virtual ~RooTsallis() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x ;
  RooRealProxy m ;
  RooRealProxy n ;
  RooRealProxy T ;
 
  Double_t evaluate() const ;

private:

 ClassDef(RooTsallis,1) // Your description goes here...
};
 
#endif
