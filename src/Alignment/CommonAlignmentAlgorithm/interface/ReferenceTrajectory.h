#ifndef REFERENCE_TRAJECTORY_H
#define REFERENCE_TRAJECTORY_H

/**
 * Author     : Gero Flucke (based on code by Edmund Widl replacing ORCA's TkReferenceTrack)
 * date       : 2006/09/17
 * last update: $Date: 2006/10/10 16:32:23 $
 * by         : $Author: ewidl $
 *
 *  Class implementing the reference trajectory of a single charged
 *  particle, i.e. a helix with 5 parameters. Given the
 *  TrajectoryStateOnSurface at the first hit and the list of all hits
 *  the local measurements, derivatives etc. as described in (and
 *  accessed via) ReferenceTrajectoryBase are calculated.
 * 
 *  The covariance-matrix of the measurements may include effects of
 *  multiple-scattering or energy-loss effects or both. This can be
 *  defined in the constructor via the variable 'materialEffects
 *  (cf. ReferenceTrajectoryBase):
 *
 *  materialEffects =  none/multipleScattering/energyLoss/combined
 *
 *  By default, the mass is assumed to be the muon-mass, but can be
 *  changed via a constructor argument.
 *
 * LIMITATIONS:
 *  So far all input hits have to be valid, but invalid hits
 *  would be needed to take into account the material effects in them...
 *
 */

#include "Alignment/CommonAlignmentAlgorithm/interface/ReferenceTrajectoryBase.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

class TrajectoryStateOnSurface;
class MagneticField;
class MaterialEffectsUpdator;
class BoundPlane;

class ReferenceTrajectory : public ReferenceTrajectoryBase
{

public:
  /**Constructor with Tsos at first hit (in physical order) and list of hits 
     [if (hitsAreReverse) ==> order of hits is in opposite direction compared
     to the flight of particle, but note that ReferenceTrajectory::recHits()
     returns the hits always in order of flight],
     the material effects to be considered and a particle mass,
     the magnetic field of the event is needed for propagations etc.
   */
  ReferenceTrajectory(const TrajectoryStateOnSurface &referenceTsos,
		      const TransientTrackingRecHit::ConstRecHitContainer &recHits,
		      bool hitsAreReverse,
		      const MagneticField *magField,
		      MaterialEffects materialEffects = combined, 
		      double mass = 0.10565836); // FIXME: ugly hard coded muon mass
  virtual ~ReferenceTrajectory() {}

  virtual ReferenceTrajectory* clone() const
    { return new ReferenceTrajectory(*this); }

private:
  /** internal method to calculate members
   */
  bool construct(const TrajectoryStateOnSurface &referenceTsos, 
		 const TransientTrackingRecHit::ConstRecHitContainer &recHits,
		 double mass, MaterialEffects materialEffects,
		 const MagneticField *magField);

  /** internal method to get apropriate updator
   */
  MaterialEffectsUpdator* createUpdator(MaterialEffects materialEffects, double mass) const;

  /** internal method to calculate jacobian
   */
  bool propagate(const BoundPlane &previousSurface, const TrajectoryStateOnSurface &previousTsos,
		 const BoundPlane &newSurface, TrajectoryStateOnSurface &newTsos, AlgebraicMatrix &newJacobian,
		 const MagneticField *magField) const;
  
  /** internal method to fill measurement and error matrix for hit iRow/2
   */
  void fillMeasurementAndError(const TransientTrackingRecHit::ConstRecHitPointer &hitPtr,
			       unsigned int iRow,
			       const TrajectoryStateOnSurface &updatedTsos);

  /** internal method to fill derivatives for hit iRow/2
   */
  void fillDerivatives(const AlgebraicMatrix &projection,
		       const AlgebraicMatrix &fullJacobian, unsigned int iRow);

  /** internal method to fill the trajectory positions for hit iRow/2
   */
  void fillTrajectoryPositions(const AlgebraicMatrix &projection, 
			       const AlgebraicVector &mixedLocalParams, 
			       unsigned int iRow);

  /** internal method to add material effects to measurments covariance matrix
   */
  void addMaterialEffectsCov(const std::vector<AlgebraicMatrix> &allJacobians, 
			     const std::vector<AlgebraicMatrix> &allProjections,
			     const std::vector<AlgebraicSymMatrix> &allCurvChanges,
			     const std::vector<AlgebraicSymMatrix> &allDeltaParaCovs);
};

#endif
