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

#ifndef ROO_MELABKG
#define ROO_MELABKG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TH2F.h"

class RooRealVar;
class RooAbsReal;
class TH2F;
class TFile;

 
class RooMELAModelBkg : public RooAbsPdf {
public:
  RooMELAModelBkg();
  RooMELAModelBkg(const char *name, const char *title,
		 RooAbsReal& _mZZ,
		 RooAbsReal& _D,
		 RooAbsReal& _a0,
		 RooAbsReal& _a1,
		 RooAbsReal& _a2,
		 RooAbsReal& _a3,
		 RooAbsReal& _a4,
		 RooAbsReal& _a5,
		 RooAbsReal& _a6,
		 RooAbsReal& _a7,
		 RooAbsReal& _a8,
		 RooAbsReal& _a9
		  );
  RooMELAModelBkg(const RooMELAModelBkg& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMELAModelBkg(*this,newname); }
  inline virtual ~RooMELAModelBkg() { }

protected:

  RooRealProxy mZZ ;
  RooRealProxy D ;
  RooRealProxy a0 ;
  RooRealProxy a1 ;
  RooRealProxy a2 ;
  RooRealProxy a3 ;
  RooRealProxy a4 ;
  RooRealProxy a5 ;
  RooRealProxy a6 ;
  RooRealProxy a7 ;
  RooRealProxy a8 ;
  RooRealProxy a9 ;

  Double_t evaluate() const ;

private:

  ClassDef(RooMELAModelBkg,1) 
};
 
#endif
