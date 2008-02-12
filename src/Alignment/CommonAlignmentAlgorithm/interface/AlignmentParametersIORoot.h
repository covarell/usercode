#ifndef Alignment_CommonAlignmentAlgorithm_AlignmentParametersIORoot_h
#define Alignment_CommonAlignmentAlgorithm_AlignmentParametersIORoot_h

/// \class AlignmentParametersIORoot
///
/// Concrete class for ROOT-based I/O of AlignmentParameters 
///
///  $Date: 2007/03/16 16:35:03 $
///  $Revision: 1.4 $
/// (last update by $Author: flucke $)

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentIORootBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParametersIO.h"

class Alignable;
class AlignmentParameters;

class AlignmentParametersIORoot : public AlignmentIORootBase, public AlignmentParametersIO
{
  friend class AlignmentIORoot;

  private:

  /// Constructor 
  AlignmentParametersIORoot(); 

  /// Write AlignmentParameters of one Alignable 
  int writeOne(Alignable* ali);

  /// Read AlignmentParameters of one Alignable 
  AlignmentParameters* readOne(Alignable* ali, int& ierr);

  /// Open IO 
  int open(const char* filename, int iteration, bool writemode)
    {return openRoot(filename,iteration,writemode);};

  /// Close IO 
  int close(void) {return closeRoot();};

  // helper functions

  /// Find entry number corresponding to Id. Returns -1 on failure.
  int findEntry(unsigned int detId,int comp);

  /// Create all branches and give names
  void createBranches(void);

  /// Set branch adresses
  void setBranchAddresses(void);

  // Alignment parameter tree 
  int theObjId, theCovRang, theCovarRang, theHieraLevel;
  unsigned int theId;
  double thePar[nParMax],theCov[nParMax*(nParMax+1)/2];

};

#endif
