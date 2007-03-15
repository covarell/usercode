
#include "Alignment/CommonAlignment/interface/AlignmentUserVariables.h"

class CSA06UserVariables : public AlignmentUserVariables {

  public:

  /** constructor */
  CSA06UserVariables(int npar) :
    jtvj(npar,0) , 
    jtve(npar,0) ,
    nhit(0) 
    //iterpar(maxiter,npar,0),
    //iterpos(maxiter,3,0),
    //iterrot(maxiter,9,0),
    //iterrpos(maxiter,3,0),
    //iterrrot(maxiter,9,0),
    //niter(0)  
  {}

  /** destructor */
  virtual ~CSA06UserVariables() {};

  /** data members */

  //static const int maxiter = 9;

  AlgebraicSymMatrix jtvj;
  AlgebraicVector jtve;
  int nhit;
  //AlgebraicMatrix iterpar;
  //AlgebraicMatrix iterpos,iterrot;
  //AlgebraicMatrix iterrpos,iterrrot;
  //int niter;

 /** clone method (copy constructor) */
  CSA06UserVariables* clone(void) const { 
    return new CSA06UserVariables(*this);
  }

};
