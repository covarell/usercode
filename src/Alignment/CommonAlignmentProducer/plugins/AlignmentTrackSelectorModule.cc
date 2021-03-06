
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "Alignment/CommonAlignmentProducer/interface/AlignmentTrackSelector.h"

// the following include is necessary to clone all track branches
// including recoTrackExtras and TrackingRecHitsOwned.
// if you remove it the code will compile, but the cloned
// tracks have only the recoTracks branch!
#include "PhysicsTools/RecoAlgos/interface/TrackSelector.h"

struct TrackConfigSelector {

  typedef std::vector<const reco::Track*> container;
  typedef container::const_iterator const_iterator;
  typedef reco::TrackCollection collection; 

  TrackConfigSelector( const edm::ParameterSet & cfg ) :
    theSelector(cfg) {}

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  size_t size() const { return selected_.size(); }

  void select( const edm::Handle<reco::TrackCollection> & c,  const edm::Event & evt, const edm::EventSetup& es) {

    all_.clear();
    selected_.clear();
    for (collection::const_iterator i = c.product()->begin(), iE = c.product()->end();
         i != iE; ++i){
      all_.push_back(& * i );
    }
    selected_=theSelector.select(all_,evt, es);
  }

private:
  container all_,selected_;
  AlignmentTrackSelector theSelector;
};

typedef ObjectSelector<TrackConfigSelector>  AlignmentTrackSelectorModule;

DEFINE_FWK_MODULE( AlignmentTrackSelectorModule );
