/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooHistPdfConv.cxx,v 1.1 2009/11/10 10:46:33 pellicci Exp $
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

//ClassImp(RooHistPdfConv);


//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooDataHist& datahist) :
  RooAbsPdf(name,title), 
  xIn(_xIn.GetName(),_xIn.GetTitle(),this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,(RooAbsReal&)RooRealConstant::value(1)),
  ssf("ssf","Sigma Scale Factor",this,(RooAbsReal&)RooRealConstant::value(1))
{  
  _histpdf = new RooDataHist(datahist);
  _variableName = "JpsictTrue";
}



//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _msSF, RooDataHist& datahist) : 
  RooAbsPdf(name,title), 
  xIn(_xIn.GetName(),_xIn.GetTitle(),this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_msSF),
  ssf("ssf","Sigma Scale Factor",this,_msSF)
{
  _histpdf = new RooDataHist(datahist);
  _variableName = "JpsictTrue";
}



//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const char *name, const char *title, RooAbsReal& _xIn, 
			     RooAbsReal& _mean, RooAbsReal& _sigma, 
			     RooAbsReal& _meanSF, RooAbsReal& _sigmaSF,
                             RooDataHist& datahist ) : 
  RooAbsPdf(name,title), 
  xIn(_xIn.GetName(),_xIn.GetTitle(),this,_xIn),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma),
  msf("msf","Mean Scale Factor",this,_meanSF),
  ssf("ssf","Sigma Scale Factor",this,_sigmaSF)
{   
  _histpdf = new RooDataHist(datahist); 
  _variableName = "JpsictTrue";
}   


//_____________________________________________________________________________
RooHistPdfConv::RooHistPdfConv(const RooHistPdfConv& other, const char* name) : 
  RooAbsPdf(other,name),
  xIn(other.xIn.GetName(),this,other.xIn),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma),
  msf("msf",this,other.msf),
  ssf("ssf",this,other.ssf)
{
  _histpdf = other._histpdf;
  _variableName = other._variableName;
}


//_____________________________________________________________________________
Double_t RooHistPdfConv::evaluate() const 
{  
  // cout << "RooHistPdfConv::evaluate(" << GetName() << ")" << endl ;
  
  static const Double_t root2(sqrt(2.)) ; 
  const RooArgSet* aRow;
  RooRealVar* xprime;
 
  Double_t result(0) ;
  // *** Convolution with hist PDF ***
  for (Int_t i=0; i<_histpdf->numEntries(); i++) {
    
    aRow = _histpdf->get(i);
    Double_t halfBinSize = _histpdf->binVolume(*aRow)/2.0;
    Double_t weight = _histpdf->weight(*aRow,0,false)/_histpdf->sum(false);
    xprime = (RooRealVar*)aRow->find(_variableName.c_str());

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
  const RooArgSet* aRow;
  RooRealVar* xprime;

  Double_t result(0.);
  for (Int_t i=0; i<_histpdf->numEntries(); i++) {
    aRow = _histpdf->get(i);
    Double_t halfBinSize = _histpdf->binVolume(*aRow)/2.0;
    Double_t weight = _histpdf->weight(*aRow,0,false)/_histpdf->sum(false);
    xprime = (RooRealVar*)aRow->find(_variableName.c_str());

    result += 0.5*weight*(cerfIndefiniteInt(xprime->getVal() - halfBinSize) - cerfIndefiniteInt(xprime->getVal() + halfBinSize) );
  }

  return result;
}

Double_t RooHistPdfConv::cerfIndefiniteInt(Double_t xi) const
{
  static const Double_t root2(sqrt(2.));
  const Double_t a = -1./(root2*sigma*ssf);
  const Double_t b = xi + mean*ssf;

  return -b/a*RooMath::erf(b+a*xIn) + xIn*RooMath::erfc(b+a*xIn) - (exp(-b*b-2*a*b*xIn-a*a*xIn*xIn));
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




