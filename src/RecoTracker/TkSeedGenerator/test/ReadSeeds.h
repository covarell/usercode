#ifndef RecoTracker_TkSeedGenerator_ReadSeeds_h
#define RecoTracker_TkSeedGenerator_ReadSeeds_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


class ReadSeeds : public edm::EDAnalyzer
{
 public:
  
  explicit ReadSeeds(const edm::ParameterSet& conf);
  
  virtual ~ReadSeeds();
  
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
 private:
  edm::ParameterSet conf_;
  
};


#endif
