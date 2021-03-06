#ifndef Alignment_CommonAlignmentAlgorithm_AlignmentCorrelationsStore_h
#define Alignment_CommonAlignmentAlgorithm_AlignmentCorrelationsStore_h

#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"

#include <map>

/// This class manages the storage and retrieval of correlations between Alignables
/// for the AlignmentParameterStore. This basic implementation simply stores the
/// offdiagonal entries of the "big covariance matrix". ATTENTION: Definition of
/// data structure "Correlations" differs from definition in AlignmentParameterStore,
/// but is used only internally.

class AlignmentCorrelationsStore
{

public:

  typedef std::map< Alignable*, AlgebraicMatrix > CorrelationsTable;
  typedef std::map< Alignable*, CorrelationsTable* > Correlations;

  AlignmentCorrelationsStore( void );

  virtual ~AlignmentCorrelationsStore( void ) {}

  /// Write correlations directly to the covariance matrix starting at the
  /// given position. Indices are assumed to start from 0.
  virtual void correlations( Alignable* ap1, Alignable* ap2,
			     AlgebraicSymMatrix& cov, int row, int col ) const;

  /// Get correlations directly from the given position of the covariance
  /// matrix and store them. Indices are assumed to start from 0.
  virtual void setCorrelations( Alignable* ap1, Alignable* ap2,
				const AlgebraicSymMatrix& cov, int row, int col );

  /// Set correlations.
  virtual void setCorrelations( Alignable* ap1, Alignable* ap2,	AlgebraicMatrix& mat );

  /// Check whether correlations are stored for a given pair of alignables.
  virtual bool correlationsAvailable( Alignable* ap1, Alignable* ap2 ) const;

  /// Reset correlations.
  virtual void resetCorrelations( void );

  /// Get number of stored correlations.
  virtual unsigned int size( void ) const;

protected:

  virtual void fillCorrelationsTable( Alignable* ap1, Alignable* ap2, CorrelationsTable* table,
				      const AlgebraicSymMatrix& cov, int row, int col, bool transpose );

  virtual void fillCovariance( Alignable* ap1, Alignable* ap2, const AlgebraicMatrix& entry,
			       AlgebraicSymMatrix& cov, int row, int col ) const;

  virtual void fillCovarianceT( Alignable* ap1, Alignable* ap2, const AlgebraicMatrix& entry,
				AlgebraicSymMatrix& cov, int row, int col ) const;

  virtual void readFromCovariance( Alignable* ap1, Alignable* ap2, AlgebraicMatrix& entry,
				   const AlgebraicSymMatrix& cov, int row, int col );

  virtual void readFromCovarianceT( Alignable* ap1, Alignable* ap2, AlgebraicMatrix& entry,
				    const AlgebraicSymMatrix& cov, int row, int col );

  Correlations theCorrelations;

};

#endif
