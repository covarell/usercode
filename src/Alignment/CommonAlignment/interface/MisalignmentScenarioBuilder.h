#ifndef Alignment_CommonAlignment_MisalignmentScenarioBuilder_h
#define Alignment_CommonAlignment_MisalignmentScenarioBuilder_h

/// \class MisalignmentScenarioBuilder
///
/// $Date: 2007/01/12 09:47:39 $
/// $Revision: 1.1 $
///
/// $Author: fronga $
/// \author Frederic Ronga - CERN-PH-CMG

#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableModifier.h"
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"

/// Base class to build a scenario from configuration and apply to either tracker or muon.

class MisalignmentScenarioBuilder
{

public:
 
  /// Default constructor
  MisalignmentScenarioBuilder( ) {};

  /// Destructor
  virtual ~MisalignmentScenarioBuilder() {};

  /// Apply misalignment scenario to the tracker (sub-system specific)
  virtual void applyScenario( const edm::ParameterSet& scenario ) = 0;

protected: // Methods

  /// Decode movements defined in given parameter set for given set of alignables
  void decodeMovements_( const edm::ParameterSet& pSet, std::vector<Alignable*> alignables );
  
  /// Decode movements defined in given parameter set for given set of alignables tagged by given name
  void decodeMovements_( const edm::ParameterSet& pSet, std::vector<Alignable*> alignables,
						 std::string levelName );

  /// Apply movements given by parameter set to given alignable
  void applyMovements_( Alignable* alignable, const edm::ParameterSet& pSet );
  
  /// Merge two sets of parameters into one (the first argument)
  void mergeParameters_( edm::ParameterSet& localSet, const edm::ParameterSet& globalSet ) const;

  /// Propagate global parameters to sub-parameters
  void propagateParameters_( const edm::ParameterSet& pSet, const std::string& globalName,
							 edm::ParameterSet& subSet ) const;

  /// Get parameter set corresponding to given name (returns empty parameter set if does not exist)
  edm::ParameterSet getParameterSet_( const std::string& name, const edm::ParameterSet& pSet ) const;

  /// Check if given parameter exists in parameter set
  bool hasParameter_( const std::string& name, const edm::ParameterSet& pSet ) const;

  /// Print all parameters and values for given set
  void printParameters_( const edm::ParameterSet& pSet, const bool showPsets = false ) const;

  /// Check if given parameter is for a top-level structure
  const bool isTopLevel_( const std::string& parameterSetName ) const; 

  /// Get root name of a parameter set (e.g. 'Rod' in 'Rods' or 'Rod1')
  const std::string rootName_( const std::string& parameterSetName ) const;
  

protected: // Members

  edm::ParameterSet theScenario;           ///< Misalignment scenario to apply (from config file)
  AlignableModifier theModifier;           ///< Helper class for random movements
  
  AlignableObjectId theAlignableObjectId;  ///< Type to name converter
  
  int theModifierCounter;                  ///< Counter for applied modification

  std::string indent;                      ///< Depth in hierarchy
  

};



#endif
