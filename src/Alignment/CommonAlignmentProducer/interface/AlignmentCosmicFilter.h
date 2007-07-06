
#ifndef ALIGNMENTCOSMICFILTER_H
#define ALIGNMENTCOSMICFILTER_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "TFile.h"
#include "TH1F.h"

#include <map>
#include <vector>

using namespace std;

class AlignmentCosmicFilter : public edm::EDFilter {
  public:
  explicit AlignmentCosmicFilter(const edm::ParameterSet& conf);
  virtual ~AlignmentCosmicFilter() {}
  //   virtual bool filter(edm::Event & e, edm::EventSetup const& c);
  virtual void beginJob(edm::EventSetup const& c);
  virtual void endJob();  
  bool filter(edm::Event & iEvent, edm::EventSetup const& c);
  

 private:
  edm::ParameterSet conf_;

  bool applySeedNumber;
  int minNSeeds,maxNSeeds;

  bool applyBasicCuts,applyMultiplicityFilter;
  int minMultiplicity,maxMultiplicity;
  double ptMin,ptMax,etaMin,etaMax,phiMin,phiMax,nHitMin,nHitMax,chi2nMax;
  int minHitsinTIB, minHitsinTOB, minHitsinTID, minHitsinTEC;
  
  TrackerAlignableId *TkMap;

  // test
  TFile* testFile;
  TH1F* hnSeeds;
  TH1F* hnTracks;
  TH1F* hPt;
  TH1F* hPhi;
  TH1F* hEta;
  TH1F* hnHits;
  TH1F* hChi2;
  TH1F* hnHitsinTIB;
  TH1F* hnHitsinTOB;
  TH1F* hnHitsinTID;
  TH1F* hnHitsinTEC;

};

#endif 
