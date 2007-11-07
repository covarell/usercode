#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Alignment/CommonAlignmentProducer/interface/AlignmentTrackSelector.h"

// constructor ----------------------------------------------------------------

AlignmentTrackSelector::AlignmentTrackSelector(const edm::ParameterSet & cfg) :
  conf_(cfg),
  applyBasicCuts( cfg.getParameter<bool>( "applyBasicCuts" ) ),
  applyNHighestPt( cfg.getParameter<bool>( "applyNHighestPt" ) ),
  applyMultiplicityFilter( cfg.getParameter<bool>( "applyMultiplicityFilter" ) ),
  seedOnlyFromAbove( cfg.getParameter<bool>( "seedOnlyFromAbove" ) ),
  applyIsolation( cfg.getParameter<bool>( "applyIsolationCut" ) ),
  chargeCheck( cfg.getParameter<bool>( "applyChargeCheck" ) ),
  nHighestPt( cfg.getParameter<int>( "nHighestPt" ) ),
  minMultiplicity ( cfg.getParameter<int>( "minMultiplicity" ) ),
  maxMultiplicity ( cfg.getParameter<int>( "maxMultiplicity" ) ),
  ptMin( cfg.getParameter<double>( "ptMin" ) ),
  ptMax( cfg.getParameter<double>( "ptMax" ) ),
  etaMin( cfg.getParameter<double>( "etaMin" ) ),
  etaMax( cfg.getParameter<double>( "etaMax" ) ),
  phiMin( cfg.getParameter<double>( "phiMin" ) ),
  phiMax( cfg.getParameter<double>( "phiMax" ) ),
  nHitMin( cfg.getParameter<double>( "nHitMin" ) ),
  nHitMax( cfg.getParameter<double>( "nHitMax" ) ),
  chi2nMax( cfg.getParameter<double>( "chi2nMax" ) ),
  isoCut( cfg.getParameter<double>( "isolationCut" ) ),
  chargeCut( cfg.getParameter<double>( "chargeCut" ) )
  
{

  if (applyBasicCuts)
	edm::LogInfo("AlignmentTrackSelector") 
	  << "applying basic track cuts ..."
	  << "\nptmin,ptmax:     " << ptMin   << "," << ptMax 
	  << "\netamin,etamax:   " << etaMin  << "," << etaMax
	  << "\nphimin,phimax:   " << phiMin  << "," << phiMax
	  << "\nnhitmin,nhitmax: " << nHitMin << "," << nHitMax
	  << "\nchi2nmax:        " << chi2nMax;

  if (applyNHighestPt)
	edm::LogInfo("AlignmentTrackSelector") 
	  << "filter N tracks with highest Pt N=" << nHighestPt;

  if (applyMultiplicityFilter)
	edm::LogInfo("AlignmentTrackSelector") 
	  << "apply multiplicity filter N >= " << minMultiplicity << " and N <= " << maxMultiplicity;

  edm::ParameterSet minHitsPerSubdet = conf_.getParameter<edm::ParameterSet>( "minHitsPerSubDet" );
  minHitsinTIB = minHitsPerSubdet.getParameter<int>( "inTIB" );
  minHitsinTOB = minHitsPerSubdet.getParameter<int>( "inTOB" );
  minHitsinTID = minHitsPerSubdet.getParameter<int>( "inTID" );
  minHitsinTEC = minHitsPerSubdet.getParameter<int>( "inTEC" );
  
  edm::LogInfo("AlignmentTrackSelector") 
    << "Minimum number of hits in TIB/TID/TOB/TEC = " << minHitsinTIB << "/" << minHitsinTID << "/" << minHitsinTOB << "/" << minHitsinTEC;

  TkMap = new TrackerAlignableId();
}

// destructor -----------------------------------------------------------------

AlignmentTrackSelector::~AlignmentTrackSelector()
{}


// do selection ---------------------------------------------------------------

AlignmentTrackSelector::Tracks 
AlignmentTrackSelector::select(const Tracks& tracks, const edm::Event& evt, const edm::EventSetup& es ) const 
{
  Tracks result=tracks;

  // apply basic track cuts (if selected)
  if (applyBasicCuts)  result= this->basicCuts(result, evt, es);

  // filter N tracks with highest Pt (if selected)
  if (applyNHighestPt) result= this->theNHighestPtTracks(result);

  // apply minimum multiplicity requirement (if selected)
  if (applyMultiplicityFilter) {
    if (result.size()<(unsigned int)minMultiplicity || result.size()>(unsigned int)maxMultiplicity ) result.clear();
  }

  //edm::LogDebug("AlignmentTrackSelector") << "tracks all,kept: " << tracks.size() << "," << result.size();

  return result;

}

// make basic cuts ------------------------------------------------------------

AlignmentTrackSelector::Tracks 
AlignmentTrackSelector::basicCuts(const Tracks& tracks, const edm::Event& evt, const edm::EventSetup& es ) const 
{
  Tracks result;

  for(Tracks::const_iterator it=tracks.begin();
      it!=tracks.end();it++) {
    const reco::Track* trackp=*it;
    float pt=trackp->pt();
    float eta=trackp->eta();
    float phi=trackp->phi();
    int nhit = trackp->numberOfValidHits(); 
    float chi2n = trackp->normalizedChi2();
    int nhitinTIB = 0;
    int nhitinTOB = 0;
    int nhitinTID = 0;
    int nhitinTEC = 0;

    //edm::LogDebug("AlignmentTrackSelector") << " pt,eta,phi,nhit: "
    //  <<pt<<","<<eta<<","<<phi<<","<<nhit;

    if (pt>ptMin && pt<ptMax 
       && eta>etaMin && eta<etaMax 
       && phi>phiMin && phi<phiMax 
       && nhit>=nHitMin && nhit<=nHitMax
       && chi2n<chi2nMax) {
 
      int thishit = 0;
      bool okSeed = true;
      bool okIso = true; 
      bool okCharge = true;

      for (trackingRecHit_iterator iHit = trackp->recHitsBegin(); iHit != trackp->recHitsEnd(); iHit++) {
        thishit++;
	std::pair<int,int> typeAndLay = TkMap->typeAndLayerFromDetId( (*iHit)->geographicalId() ); 
	int type = typeAndLay.first; 

        if (seedOnlyFromAbove && thishit == 1 && type == int(StripSubdetector::TOB)) okSeed = false;  
	// if first hit is in TOB seed is from below (mysteries of tracking...)
           
        if ((*iHit)->isValid()) {

	  if (chargeCheck) {
            const TrackingRecHit* therechit = (*iHit).get();
            float charge1 = 0;
            float charge2 = 0;
	    const SiStripMatchedRecHit2D* matchedhit = dynamic_cast<const SiStripMatchedRecHit2D*>(therechit);
	    const SiStripRecHit2D* hit = dynamic_cast<const SiStripRecHit2D*>(therechit);
	    
	    if (matchedhit) {  
	      const SiStripRecHit2D *monohit=matchedhit->monoHit();    
	      const SiStripCluster* monocluster = &*(monohit->cluster());
	      const std::vector<uint16_t> amplitudesmono( monocluster->amplitudes().begin(),
							  monocluster->amplitudes().end());
	      for(size_t ia=0; ia<amplitudesmono.size();ia++)
		{ charge1+=amplitudesmono[ia];} 
	      
              const SiStripRecHit2D *stereohit=matchedhit->stereoHit();   
	      const SiStripCluster* stereocluster = &*(stereohit->cluster());
	      const std::vector<uint16_t> amplitudesstereo( stereocluster->amplitudes().begin(),
							    stereocluster->amplitudes().end());
	      for(size_t ia=0; ia<amplitudesstereo.size();ia++)
		{charge2+=amplitudesstereo[ia];}
	      // std::cout << "charge1 = " << charge1 << "\n";
	      // std::cout << "charge2 = " << charge2 << "\n";
	      if (charge1 < chargeCut || charge2 < chargeCut) okCharge = false;
	    }
	    else if (hit) {
	      
	      const SiStripCluster* cluster = &*(hit->cluster());
	      const std::vector<uint16_t> amplitudes( cluster->amplitudes().begin(),
						      cluster->amplitudes().end());
	      for(size_t ia=0; ia<amplitudes.size();ia++)
		{charge1+=amplitudes[ia];}
	      // std::cout << "charge1 = " << charge1 << "\n";
	      if (charge1 < chargeCut) okCharge = false;
	    }
	  }
          
	  if (applyIsolation) {

            edm::ESHandle<TrackerGeometry> tracker;
	    es.get<TrackerDigiGeometryRecord>().get(tracker);

	    edm::Handle<SiStripRecHit2DCollection> rphirecHits;
	    edm::InputTag rphirecHitsTag = conf_.getParameter<edm::InputTag>("rphirecHits");
	    evt.getByLabel( rphirecHitsTag, rphirecHits );
	    
	    edm::Handle<SiStripMatchedRecHit2DCollection> matchedrecHits;
	    edm::InputTag matchedrecHitsTag = conf_.getParameter<edm::InputTag>("matchedrecHits");
	    evt.getByLabel( matchedrecHitsTag, matchedrecHits );
	    
	    SiStripRecHit2DCollection::const_iterator istripSt; 
	    SiStripMatchedRecHit2DCollection::const_iterator istripStm; 
	    SiStripRecHit2DCollection stripcollSt = *rphirecHits;
	    SiStripMatchedRecHit2DCollection stripcollStm = *matchedrecHits;
	    
            DetId idet = (*iHit)->geographicalId(); 
            GlobalPoint point = tracker->idToDet(idet)->surface().toGlobal((*iHit)->localPosition());
 
	    for( istripSt=stripcollSt.begin(); istripSt!=stripcollSt.end(); istripSt++ ) {
	      const SiStripRecHit2D *aHit = &*(istripSt);
              DetId mydet1 = aHit->geographicalId(); 
	      GlobalPoint point1 = tracker->idToDet(mydet1)->surface().toGlobal(aHit->localPosition());
	      float theDistance = sqrt( pow(point1.x() - point.x(), 2) +
					pow(point1.y() - point.y(), 2) +
					pow(point1.z() - point.z(), 2) );
	      // std::cout << "theDistance1 = " << theDistance << "\n";
	      if (idet.rawId() == mydet1.rawId() && theDistance > 0.001 && theDistance < isoCut) okIso = false;
	    }
	    
	    for( istripStm=stripcollStm.begin(); istripStm!=stripcollStm.end(); istripStm++ ) {
	      const SiStripMatchedRecHit2D *aHit = &*(istripStm);
	      DetId mydet2 = aHit->geographicalId(); 
	      GlobalPoint point2 = tracker->idToDet(mydet2)->surface().toGlobal(aHit->localPosition());
	      float theDistance = sqrt( pow(point2.x() - point.x(), 2) +
					pow(point2.y() - point.y(), 2) +
					pow(point2.z() - point.z(), 2) );
              // std::cout << "theDistance1 = " << theDistance << "\n";
	      if (idet.rawId() == mydet2.rawId() && theDistance > 0.001 && theDistance < isoCut) okIso = false;
	    }
	  }     
            
	  
	  if (type == int(StripSubdetector::TIB)) nhitinTIB++;
	  if (type == int(StripSubdetector::TOB)) nhitinTOB++;
	  if (type == int(StripSubdetector::TID)) nhitinTID++;
	  if (type == int(StripSubdetector::TEC)) nhitinTEC++;
	}
      }
      
      if (nhitinTIB>=minHitsinTIB &&
	  nhitinTOB>=minHitsinTOB &&
	  nhitinTID>=minHitsinTID &&
	  nhitinTEC>=minHitsinTEC &&
          okSeed && okIso && okCharge) result.push_back(trackp);
    }
  }
  
  return result;
}

//-----------------------------------------------------------------------------

AlignmentTrackSelector::Tracks 
AlignmentTrackSelector::theNHighestPtTracks(const Tracks& tracks) const
{
  Tracks sortedTracks=tracks;
  Tracks result;

  // sort in pt
  std::sort(sortedTracks.begin(),sortedTracks.end(),ptComparator);

  // copy theTrackMult highest pt tracks to result vector
  int n=0;
  for (Tracks::const_iterator it=sortedTracks.begin();
	   it!=sortedTracks.end(); it++) {
	if (n<nHighestPt) { result.push_back(*it); n++; }
  }

  return result;
}

