#ifndef CosmicTrackFinder_h
#define CosmicTrackFinder_h

// Package:    RecoTracker/SingleTrackPattern
// Class:      CosmicTrackFinder
// Original Author:  Michele Pioppi-INFN perugia


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "RecoTracker/SingleTrackPattern/interface/CosmicTrajectoryBuilder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

namespace cms
{
  class CompareTrajLay {
  public:
    bool operator()(Trajectory *t1,
		    Trajectory *t2){
      AnalHits(t1->recHits());
      uint alay=nlay;
      AnalHits(t2->recHits());
      uint blay=nlay;
      if (alay!=blay) return alay > blay;
      if (t1->foundHits() != t2->foundHits()) 
	return t1->foundHits()> t2->foundHits();
      return t1->chiSquared()< t2->chiSquared();
      // std::cout<<"chi "<<t1.chiSquared()<<" "<<t2.chiSquared()<<std::endl;
      // return false;
    }
    void  AnalHits(std::vector< ConstReferenceCountingPointer< TransientTrackingRecHit> > hits){
      ltob1=false; ltob2=false; ltib1=false; ltib2=false;
      std::vector< ConstReferenceCountingPointer< TransientTrackingRecHit> >::const_iterator hit;
      //     ConstRecHitIterator hit;
      for(hit=hits.begin();hit!=hits.end();hit++){
	uint iid=(*hit)->hit()->geographicalId().rawId();
	
	int sub=(iid>>25)&0x7 ;
	int lay=(iid>>16) & 0xF;
	if ((lay==1)&&(sub==3)) ltib1=true;
	if ((lay==2)&&(sub==3)) ltib2=true;
	if ((lay==1)&&(sub==5)) ltob1=true;
	if ((lay==2)&&(sub==5)) ltob2=true;
    }
      nlay=ltib1+ltib2+ltob1+ltob2;
      
    }
    
  private:
    bool ltib1,ltib2,ltob1,ltob2;
    uint nlay;
    
  };
  class CompareTrajChi {
  public:
    bool operator()(Trajectory *t1,
		    Trajectory *t2){
      if (t1->foundHits() != t2->foundHits()) 
	return t1->foundHits()> t2->foundHits();
      return t1->chiSquared()< t2->chiSquared();  
    }
  };
  class CosmicTrackFinder : public edm::EDProducer
  {

    typedef TrajectoryStateOnSurface     TSOS;
  public:

    explicit CosmicTrackFinder(const edm::ParameterSet& conf);

    virtual ~CosmicTrackFinder();

    virtual void produce(edm::Event& e, const edm::EventSetup& c);

  private:
    CosmicTrajectoryBuilder cosmicTrajectoryBuilder_;
    edm::ParameterSet conf_;
    std::string geometry;
    bool trinevents;
  };
}

#endif
