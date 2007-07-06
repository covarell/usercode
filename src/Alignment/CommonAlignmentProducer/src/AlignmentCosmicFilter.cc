
#include "Alignment/CommonAlignmentProducer/interface/AlignmentCosmicFilter.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include <map>
#include <vector>

using namespace std;

AlignmentCosmicFilter::AlignmentCosmicFilter(const edm::ParameterSet& conf):      conf_(conf),
  applySeedNumber( conf_.getParameter<bool>( "applySeedNumber" ) ),
  minNSeeds ( conf_.getParameter<int>( "minNSeeds" ) ),
  maxNSeeds ( conf_.getParameter<int>( "maxNSeeds" ) ),
  applyBasicCuts( conf_.getParameter<bool>( "applyBasicCuts" ) ),
  applyMultiplicityFilter( conf_.getParameter<bool>( "applyMultiplicityFilter" ) ),
  minMultiplicity ( conf_.getParameter<int>( "minMultiplicity" ) ),
  maxMultiplicity ( conf_.getParameter<int>( "maxMultiplicity" ) ),
  ptMin( conf_.getParameter<double>( "ptMin" ) ),
  ptMax( conf_.getParameter<double>( "ptMax" ) ),
  etaMin( conf_.getParameter<double>( "etaMin" ) ),
  etaMax( conf_.getParameter<double>( "etaMax" ) ),
  phiMin( conf_.getParameter<double>( "phiMin" ) ),
  phiMax( conf_.getParameter<double>( "phiMax" ) ),
  nHitMin( conf_.getParameter<double>( "nHitMin" ) ),
  nHitMax( conf_.getParameter<double>( "nHitMax" ) ),
  chi2nMax( conf_.getParameter<double>( "chi2nMax" ) )
{

  if (applySeedNumber)
	std::cout 
	  << "apply seedNumber N<=" << minNSeeds << " and N>= " << maxNSeeds 
          << std::endl;

  if (applyBasicCuts)
        std::cout 
	  << "applying basic track cuts ..."
	  << "\nptmin,ptmax:     " << ptMin   << "," << ptMax 
	  << "\netamin,etamax:   " << etaMin  << "," << etaMax
	  << "\nphimin,phimax:   " << phiMin  << "," << phiMax
	  << "\nnhitmin,nhitmax: " << nHitMin << "," << nHitMax
	  << "\nchi2nmax:        " << chi2nMax << std::endl;

  if (applyMultiplicityFilter)
	std::cout 
	  << "apply multiplicity filter N>=" << minMultiplicity << " and N>= " << maxMultiplicity
        << std::endl;

  edm::ParameterSet minHitsPerSubdet = conf_.getParameter<edm::ParameterSet>( "minHitsPerSubDet" );
  minHitsinTIB = minHitsPerSubdet.getUntrackedParameter<int>( "inTIB" , 0 );
  minHitsinTOB = minHitsPerSubdet.getUntrackedParameter<int>( "inTOB" , 0 );
  minHitsinTID = minHitsPerSubdet.getUntrackedParameter<int>( "inTID" , 0 );
  minHitsinTEC = minHitsPerSubdet.getUntrackedParameter<int>( "inTEC" , 0 );

  std::cout 
    << "Minimum # of hits in TIB/TOB/TID/TEC" << minHitsinTIB << "/" << minHitsinTOB << "/" << minHitsinTID << "/" << minHitsinTEC 
    << std::endl; 
}

void AlignmentCosmicFilter::beginJob(const edm::EventSetup& iSetup) {

  // test Rootfile
  testFile   = new TFile( "test.root", "RECREATE" ) ; 
  hnSeeds = new TH1F( "hnSeeds", "Number of seeds", 100,  0., 500. ) ;
  hnTracks = new TH1F( "hnTracks", "Number of tracks", 3,  0., 2. ) ;
  hPt = new TH1F( "hPt", "Pt track", 100,  0., 30. ) ;
  hPhi = new TH1F( "hPhi", "Phi track", 100, -3.14, 0. ) ;
  hEta = new TH1F( "hEta", "Eta track", 100,  -3.0, 3.0 ) ;
  hnHits = new TH1F( "hnHits", "Number of hits", 100,  0., 30. ) ;
  hChi2 = new TH1F( "hChi2", "Chi2 track", 100,  0., 2000. ) ;
  hnHitsinTIB = new TH1F( "hnHitsinTIB", "Number of hits in TIB", 100,  0., 20. ) ;
  hnHitsinTOB = new TH1F( "hnHitsinTOB", "Number of hits in TOB", 100,  0., 20. ) ;
  hnHitsinTID = new TH1F( "hnHitsinTID", "Number of hits in TID", 100,  0., 20. ) ;
  hnHitsinTEC = new TH1F( "hnHitsinTEC", "Number of hits in TEC", 100,  0., 20. ) ;

}

bool AlignmentCosmicFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<TrajectorySeedCollection> theSeeds;
  edm::InputTag seedTag = conf_.getParameter<edm::InputTag>("seedTag");
  iEvent.getByLabel( seedTag, theSeeds );

  edm::Handle<reco::TrackCollection> theTracks;
  edm::InputTag trackTag = conf_.getParameter<edm::InputTag>("trackTag");
  iEvent.getByLabel( trackTag, theTracks );

  hnSeeds->Fill( theSeeds.product()->size() );
  hnTracks->Fill( theTracks.product()->size() );

  if (applySeedNumber) {
    if (theSeeds.product()->size()<(unsigned int)minNSeeds || theSeeds.product()->size()>(unsigned int)maxNSeeds ) return false;
  }
  
  if (applyMultiplicityFilter) {
    if (theTracks.product()->size()<(unsigned int)minMultiplicity || theTracks.product()->size()>(unsigned int)maxMultiplicity ) return false;
  }

  if (applyBasicCuts) {

    bool thisTrackOk = false;
    for (reco::TrackCollection::const_iterator it=theTracks.product()->begin();
	it!=theTracks.product()->end();it++) {

      if (!thisTrackOk) {
	const reco::Track* trackp= &(*it);
	float pt=trackp->pt();
	float eta=trackp->eta();
	float phi=trackp->phi();
	int nhit = trackp->numberOfValidHits(); 
	float chi2n = trackp->normalizedChi2();

	hPt->Fill( pt );
	hPhi->Fill( phi );
	hEta->Fill( eta );
        hnHits->Fill( nhit );
        hChi2->Fill( chi2n );

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
	  
	  for (trackingRecHit_iterator iHit = trackp->recHitsBegin(); iHit != trackp->recHitsEnd(); iHit++) {
	    std::pair<int,int> typeAndLay = TkMap->typeAndLayerFromDetId( (*iHit)->geographicalId() );
	    int type = typeAndLay.first; 
	    if (type == int(StripSubdetector::TIB)) nhitinTIB++;
	    if (type == int(StripSubdetector::TOB)) nhitinTOB++;
	    if (type == int(StripSubdetector::TID)) nhitinTID++;
	    if (type == int(StripSubdetector::TEC)) nhitinTEC++;
	  }
	  
          hnHitsinTIB->Fill( nhitinTIB );
          hnHitsinTOB->Fill( nhitinTOB );
          hnHitsinTID->Fill( nhitinTID );
          hnHitsinTEC->Fill( nhitinTEC );
 
	  if (nhitinTIB>=minHitsinTIB &&
	      nhitinTOB>=minHitsinTOB &&
	      nhitinTID>=minHitsinTID &&
	      nhitinTEC>=minHitsinTEC ) thisTrackOk = true;
	}
      }
    }
    return thisTrackOk;
  }

  return true;
}

void AlignmentCosmicFilter::endJob() {
  testFile->Write();
  testFile->Close();
}
