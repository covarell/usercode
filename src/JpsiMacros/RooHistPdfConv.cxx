/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooHistPdfConv.cxx,v 1.9 2011/01/07 14:36:29 covarell Exp $
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

#include "TMath.h"
#include <vector>

#include "RooFit.h"
#include "Riostream.h"
#include "RooHistPdfConv.h"
#include "RooRealConstant.h"
#include "RooRandom.h"
#include "RooRealVar.h"

//ClassImp(RooHistPdfConv);

static const Double_t root2(sqrt(2.));
static const Double_t rootpi(sqrt(acos(-1.)));
static unsigned int nbins;
static std::vector<Double_t> halfBinWidths;
static std::vector<Double_t> weights;
static std::vector<Double_t> centers;

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
  _variableName = "Jpsi_CtTrue";
  init();
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
  _variableName = "Jpsi_CtTrue";
  init();
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
  _variableName = "Jpsi_CtTrue";
  init();
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
  init();
}

//_____________________________________________________________________________
void RooHistPdfConv::init() const 
{ 
  const RooArgSet* aRow;
  RooRealVar* xprime;
 
  // *** Build vectors for speed reasons ***
  nbins = _histpdf->numEntries();

  for (Int_t i=0; i<_histpdf->numEntries(); i++) {
    
    aRow = _histpdf->get(i);
    xprime = (RooRealVar*)aRow->find(_variableName.c_str());
  
    const Double_t halfBinSize = xprime->getBinning().binWidth(i)/2.0;

    halfBinWidths.push_back(halfBinSize);

    centers.push_back(xprime->getVal());

    // remove non-living components
    if ( xprime->getBinning().binLow(i)*xprime->getBinning().binHigh(i) < 0) {
      weights.push_back(0.);
    } else {
      Double_t weight = (_histpdf->weight(*aRow,0,false)/_histpdf->sum(false))*((xprime->getBinning().highBound() - xprime->getBinning().lowBound())/halfBinSize);
      weights.push_back(weight);
    }
  }
  
  return;

}

//_____________________________________________________________________________
Double_t RooHistPdfConv::evaluate() const 
{  
  // cout << "RooHistPdfConv::evaluate(" << GetName() << ")" << endl ;
 
  Double_t result(0) ;
  // *** Convolution with hist PDF ***
  for (unsigned int i=0; i<nbins; i++) {
  
    Double_t halfBinSize = halfBinWidths[i];
    Double_t center = centers[i];
    Double_t weight = weights[i];
    
    const Double_t c = (center - halfBinSize - xIn + (mean*msf)) / (root2*sigma*ssf);
    const Double_t d = (center + halfBinSize - xIn + (mean*msf)) / (root2*sigma*ssf);

    result += 0.5*weight*(TMath::Erfc(c)-TMath::Erfc(d));
  }

  return fabs(result);
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

  Double_t result(0.);
  for (unsigned int i=0; i<nbins; i++) {
    
    Double_t halfBinSize = halfBinWidths[i];
    Double_t center = centers[i];
    Double_t weight = weights[i];

    result += 0.5*weight*(cerfInt(center - halfBinSize) - cerfInt(center + halfBinSize) );
  }
  
  return result;
}


Double_t RooHistPdfConv::cerfInt(Double_t xi) const
{
  const Double_t xmin = xIn.min();
  const Double_t xmax = xIn.max();

  const Double_t a = -1./(root2*sigma*ssf);
  const Double_t b = (xi + mean*msf)/(root2*sigma*ssf);

  const Double_t maxInt = -b/a*TMath::Erf(b+a*xmax) + xmax*TMath::Erfc(b+a*xmax) - exp(-b*b-2*a*b*xmax-a*a*xmax*xmax)/(rootpi*a);
  const Double_t minInt = -b/a*TMath::Erf(b+a*xmin) + xmin*TMath::Erfc(b+a*xmin) - exp(-b*b-2*a*b*xmin-a*a*xmin*xmin)/(rootpi*a);

  return maxInt - minInt;
}

//_____________________________________________________________________________
/* Double_t RooHistPdfConv::analyticalIntegral(Int_t code, const char* rangeName) const 
{

  Double_t result(0.);

  const Double_t xmin = xIn.min();
  const Double_t xmax = xIn.max();
  const Double_t a = 1./(root2*sigma*ssf);

  for (unsigned int i=0; i<nbins; i++) {
    
    Double_t halfBinSize = halfBinWidths[i];
    Double_t center = centers[i];
    Double_t weight = weights[i];

    Double_t tmax1 = a*(center - halfBinSize - xmin + (mean*msf));
    Double_t tmin1 = a*(center - halfBinSize - xmax + (mean*msf));
    Double_t tmax2 = a*(center + halfBinSize - xmin + (mean*msf));
    Double_t tmin2 = a*(center + halfBinSize - xmax + (mean*msf));

    result += 0.5*weight*a*(cerfInt(tmax1) - cerfInt(tmin1) + cerfInt(tmax2) - cerfInt(tmin2));
  }
  
  return fabs(result);
}


Double_t RooHistPdfConv::cerfInt(Double_t xi) const
{
  return xi*TMath::Erfc(xi) - (1./rootpi)*exp(-xi*xi);
}
*/

//_____________________________________________________________________________
/* Int_t RooHistPdfConv::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
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
}*/




