#include "RecoTracker/TkSeedGenerator/test/ClusterFilter.h"

ClusterFilter::ClusterFilter(const edm::ParameterSet& conf): conf_(conf){
	produces< edm::DetSetVector<SiStripCluster> > ();
}

ClusterFilter::~ClusterFilter(){}

void ClusterFilter::beginJob(const edm::EventSetup& es){
}

void ClusterFilter::produce(edm::Event& e, const edm::EventSetup& c){
	//get inputs
	std::string stripClusterProducer = conf_.getParameter<std::string>("ClusterProducer");
   	edm::Handle<edm::DetSetVector<SiStripCluster> > clusterHandle;
   	e.getByLabel(stripClusterProducer, clusterHandle);
  	const edm::DetSetVector<SiStripCluster>* clusterCollection = clusterHandle.product();
	std::string cpe = conf_.getParameter<std::string>("StripCPE");
     	edm::ESHandle<StripClusterParameterEstimator> parameterestimator;
     	c.get<TrackerCPERecord>().get(cpe, parameterestimator); 
     	const StripClusterParameterEstimator &stripcpe(*parameterestimator);
	edm::ESHandle<TrackerGeometry> pDD;
     	c.get<TrackerDigiGeometryRecord>().get( pDD );
     	const TrackerGeometry &tracker(*pDD);	

	//ouptut collection
	std::vector< edm::DetSet<SiStripCluster> > vSiStripCluster;
	
	//loop on cluster to select the ones at positive y
	edm::DetSetVector<SiStripCluster>::const_iterator iDSV;
	for (iDSV = clusterCollection->begin(); iDSV != clusterCollection->end(); iDSV++){
		unsigned int id = iDSV->id;
		DetId detId(id);
		const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(detId);
		edm::DetSet<SiStripCluster>::const_iterator iDS;
		edm::DetSet<SiStripCluster> outputDetSet;
		outputDetSet.id = id;
		for (iDS = iDSV->data.begin(); iDS != iDSV->data.end(); iDS++){
			StripClusterParameterEstimator::LocalValues parameters=
				stripcpe.localParameters(*iDS,*stripdet);
			if (check(parameters.first, stripdet)) outputDetSet.data.push_back(*iDS);
		}
		if (outputDetSet.data.size()) vSiStripCluster.push_back(outputDetSet);
	}

	std::auto_ptr< edm::DetSetVector<SiStripCluster> > output(new edm::DetSetVector<SiStripCluster>(vSiStripCluster) );
	e.put(output);
}

GlobalPoint ClusterFilter::toGlobal(const LocalPoint& local, const StripGeomDetUnit* detunit){
	return detunit->surface().toGlobal(local);
}

bool ClusterFilter::check(const LocalPoint& local, const StripGeomDetUnit* detunit){
	if (detunit->specificType().isBarrel()) return toGlobal(local, detunit).y()>0;
	else return true;
}
