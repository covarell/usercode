#ifndef CSA06UserVariablesIORoot_H
#define CSA06UserVariablesIORoot_H

#include "TTree.h"

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentIORootBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentUserVariablesIO.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignmentParameters.h"

#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06AlignmentAlgorithm.h"

/** concrete class for ROOT based IO of AlignmentUserVariables */

class CSA06UserVariablesIORoot : public AlignmentIORootBase,
                                 public AlignmentUserVariablesIO
{

  public:

  typedef std::vector<Alignable*> Alignables;

  /** constructor */
  CSA06UserVariablesIORoot(); 

  /** write user variables */
  void writeCSA06UserVariables (const Alignables& alivec, 
    const char* filename, int iter, bool validCheck, int& ierr);

  /** read user variables */
  std::vector<AlignmentUserVariables*> readCSA06UserVariables 
    (const Alignables& alivec, const char* filename, int iter, int& ierr);



  private:

  /** write AlignmentParameters of one Alignable */
  int writeOne(Alignable* ali);

  /** read AlignmentParameters of one Alignable */
  AlignmentUserVariables* readOne(Alignable* ali, int& ierr);

  /** open IO */
  int open(const char* filename, int iteration, bool writemode)
    {newopen=true; return openRoot(filename,iteration,writemode);};

  /** close IO */
  int close(void) {return closeRoot();};

  // helper functions

  int findEntry(unsigned int detId,int comp);
  void createBranches(void);
  void setBranchAddresses(void);

  // data members

  static const int nparmax=6;

  /** alignment parameter tree */
  int ObjId;
  unsigned int Id;
  int Nhit,Nparj,Npare;
  double Jtvj[nparmax*(nparmax+1)/2];
  double Jtve[nparmax];

  bool newopen;
  typedef  std::map< std::pair<int,int> , int > treemaptype;
  treemaptype treemap;

};

#endif
