#ifndef Alignment_CSA06AlignmentAlgorithm_TrackLocalAngle_h
#define Alignment_CSA06AlignmentAlgorithm_TrackLocalAngle_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/GenericHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h" 
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"


class TrackLocalAngle 
{
 public:
	   
  explicit TrackLocalAngle(const TrackerGeometry * tracker);
  
  virtual ~TrackLocalAngle();
 
  std::pair<float,float> findtrackangle(const TrajectoryMeasurement& theTM);

  std::pair<float,float> findhitcharge(const TrajectoryMeasurement& theTM);

  std::pair<float,float> findhitbary(const TrajectoryMeasurement& theTM);

 private:
	
  const TrackerGeometry * tracker_;
	
};


#endif
