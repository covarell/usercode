#ifndef Alignment_CommonAlignmentAlgorithm_AlignmentTrackSelector_h
#define Alignment_CommonAlignmentAlgorithm_AlignmentTrackSelector_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/GenericHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include <vector>

namespace edm {
  class Event;
  class ParameterSet;
}

class TrackingRecHit;

class AlignmentTrackSelector
{

 public:

  typedef std::vector<const reco::Track*> Tracks; 

  /// constructor
  AlignmentTrackSelector(const edm::ParameterSet & cfg);

  /// destructor
  ~AlignmentTrackSelector();

  /// select tracks
  Tracks select(const Tracks& tracks, const edm::Event& evt, const edm::EventSetup& es) const;

 private:

  /// apply basic cuts on pt,eta,phi,nhit
  Tracks basicCuts(const Tracks& tracks, const edm::Event& evt) const;
  /// checking hit requirements beyond simple number of valid hits
  bool detailedHitsCheck(const reco::Track* track, const edm::Event& evt) const;
  bool isHit2D(const TrackingRecHit &hit) const;
  bool isOkCharge(const TrackingRecHit* therechit) const;
  bool isIsolated(const TrackingRecHit* therechit, const edm::Event& evt) const;

  /// filter the n highest pt tracks
  Tracks theNHighestPtTracks(const Tracks& tracks) const;

  /// compare two tracks in pt (used by theNHighestPtTracks)
  struct ComparePt {
    bool operator()( const reco::Track* t1, const reco::Track* t2 ) const {
      return t1->pt()> t2->pt();
    }
  };
  ComparePt ptComparator;

  const bool applyBasicCuts_, applyNHighestPt_, applyMultiplicityFilter_;
  const bool seedOnlyFromAbove_, applyIsolation_, chargeCheck_ ;
  const int nHighestPt_, minMultiplicity_, maxMultiplicity_;
  const bool multiplicityOnInput_; /// if true, cut min/maxMultiplicity on input instead of on final result
  const double ptMin_,ptMax_,etaMin_,etaMax_,phiMin_,phiMax_,nHitMin_,nHitMax_,chi2nMax_;
  const double minHitChargeStrip_, minHitIsolation_;
  edm::InputTag rphirecHitsTag_;
  edm::InputTag matchedrecHitsTag_;
  const unsigned int nHitMin2D_;
  const int minHitsinTIB_, minHitsinTOB_, minHitsinTID_, minHitsinTEC_, minHitsinBPIX_, minHitsinFPIX_;

};

#endif

