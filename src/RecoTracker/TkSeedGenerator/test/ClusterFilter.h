#ifndef ClusterFilter_h
#define ClusterFilter_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//Data Formats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
//RecoLocalTracker
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TrackerCPERecord.h"
//Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Vector/interface/LocalPoint.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"

class ClusterFilter : public edm::EDProducer
{
	public:
	ClusterFilter(const edm::ParameterSet& conf);	

	virtual ~ClusterFilter();

    	virtual void beginJob( const edm::EventSetup& );

    	virtual void produce(edm::Event& e, const edm::EventSetup& c);

	private:
	edm::ParameterSet conf_;
	bool check(const LocalPoint& local, const StripGeomDetUnit* detunit);
	GlobalPoint toGlobal(const LocalPoint& local, const StripGeomDetUnit* detunit);
};



#endif
