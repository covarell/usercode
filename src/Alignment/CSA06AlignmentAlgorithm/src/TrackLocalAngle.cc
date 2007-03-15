
#include <memory>
#include <string>
#include <iostream>
#include <TMath.h>
#include "Alignment/CSA06AlignmentAlgorithm/interface/TrackLocalAngle.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/Vector/interface/GlobalVector.h"
#include "Geometry/Vector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

using namespace std;
TrackLocalAngle::TrackLocalAngle(const TrackerGeometry* tracker)  
{  
  tracker_ = tracker;
}

// Virtual destructor needed.
TrackLocalAngle::~TrackLocalAngle() {  }  

std::pair<float,float> TrackLocalAngle::findtrackangle(const TrajectoryMeasurement& theTM)
{
  
  std::pair<float,float> monostereoang; 
  float angle1 = -9999.;
  float angle2 = -9999.;
  
  LogDebug("TrackLocalAngle::findtrackangle")<<"rechit and TSOS";

  TrajectoryStateOnSurface tsos = theTM.updatedState();
  const TransientTrackingRecHit::ConstRecHitPointer thit=theTM.recHit();
  const SiStripMatchedRecHit2D* matchedhit = dynamic_cast<const SiStripMatchedRecHit2D*>((*thit).hit());
  const SiStripRecHit2D* hit = dynamic_cast<const SiStripRecHit2D*>((*thit).hit());
  LocalVector trackdirection = tsos.localDirection();

  if (matchedhit) { //if matched hit...
    
    GluedGeomDet * gdet=(GluedGeomDet *)tracker_->idToDet(matchedhit->geographicalId());
    
    GlobalVector gtrkdir=gdet->toGlobal(trackdirection);
      
    //cluster and trackdirection on mono det
    
    // THIS IS THE POINTER TO THE MONO HIT OF A MATCHED HIT 
    const SiStripRecHit2D *monohit=matchedhit->monoHit();
    
    const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > monocluster=monohit->cluster();
    const GeomDetUnit * monodet=gdet->monoDet();
    
    LocalVector monotkdir=monodet->toLocal(gtrkdir);
   
    if(monotkdir.z() != 0){
      
      // THE LOCAL ANGLE (MONO)
      angle1 = atan(monotkdir.x()/monotkdir.z())*180/TMath::Pi();
    }
  
    //cluster and trackdirection on stereo det
    
    // THIS IS THE POINTER TO THE STEREO HIT OF A MATCHED HIT 
    const SiStripRecHit2D *stereohit=matchedhit->stereoHit();
    
    const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > stereocluster=stereohit->cluster();
    const GeomDetUnit * stereodet=gdet->stereoDet(); 
    LocalVector stereotkdir=stereodet->toLocal(gtrkdir);
    
    if(stereotkdir.z()!=0){
      
      // THE LOCAL ANGLE (STEREO)
      angle2 = atan(stereotkdir.x()/stereotkdir.z())*180/TMath::Pi();
      
    }
    
  }
  else if (hit) {
    //  hit= POINTER TO THE RECHIT
    const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > cluster=hit->cluster();
    
    if(trackdirection.z()!=0){
      
      // THE LOCAL ANGLE
      angle1 = atan(trackdirection.x()/trackdirection.z())*180/TMath::Pi();
    }    
  }

  monostereoang = make_pair(angle1, angle2); 
  return monostereoang;
}

