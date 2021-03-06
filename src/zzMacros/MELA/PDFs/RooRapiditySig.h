/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROORAPIDITY_SIG
#define ROORAPIDITY_SIG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooRapiditySig : public RooAbsPdf {
public:
	RooRapiditySig() {} ; 
	RooRapiditySig(const char *name, const char *title,RooAbsReal& _Y,RooAbsReal& _m,RooAbsReal& _sqrtS);
						
	RooRapiditySig(const RooRapiditySig& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooRapiditySig(*this,newname); }
	inline virtual ~RooRapiditySig() { }

	
protected:
	
	RooRealProxy Y ;
	RooRealProxy m ;
	RooRealProxy sqrtS ;

	
	Double_t evaluate() const ;
	

	
private:
	
	ClassDef(RooRapiditySig,1) // Your description goes here...
};

#endif
