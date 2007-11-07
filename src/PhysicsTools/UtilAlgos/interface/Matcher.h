#ifndef UtilAlgos_Matcher_h
#define UtilAlgos_Matcher_h
/* \class Matcher
 *
 * \author Luca Lista, INFN
 *
 */
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "PhysicsTools/UtilAlgos/interface/DeltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

namespace reco {
  namespace modules {
    template<typename C1, typename C2>
    class MatcherBase : public edm::EDProducer {
    public:
      MatcherBase( const edm::ParameterSet & );
      ~MatcherBase();

    protected:
      typedef typename C1::value_type T1;
      typedef typename C2::value_type T2;
      typedef edm::AssociationMap<edm::OneToOne<C1, C2> >MatchMap;

    private:
      void produce( edm::Event&, const edm::EventSetup& );
      edm::InputTag src_;
      edm::InputTag matched_;
      double distMin_;
      virtual double matchDistance( const T1 &, const T2 & ) const = 0;
      virtual bool select( const T1 &, const T2 & ) const = 0;
    };

    template<typename C1, typename C2,
	     typename S, typename D = DeltaR<typename C1::value_type, typename C2::value_type> >
    class Matcher : public MatcherBase<C1, C2> {
    public:
      Matcher(  const edm::ParameterSet & cfg ) : 
	MatcherBase<C1, C2>( cfg ),
        select_( reco::modules::make<S>( cfg ) ), 
	distance_( reco::modules::make<D>( cfg ) ) { }
      ~Matcher() { }
    private:
      typedef typename MatcherBase<C1, C2>::T1 T1;
      typedef typename MatcherBase<C1, C2>::T2 T2;
      typedef typename MatcherBase<C1, C2>::MatchMap MatchMap;

      double matchDistance( const T1 & c1, const T2 & c2 ) const {
	return distance_( c1, c2 );
      }
      bool select( const T1 & c1, const T2 & c2 ) const { 
	return select_( c1, c2 ); 
      }
      S select_;
      D distance_;
    };

    namespace helper {
      typedef std::pair<size_t, double> MatchPair;

      struct SortBySecond {
	bool operator()( const MatchPair & p1, const MatchPair & p2 ) const {
	  return p1.second < p2.second;
	} 
      };
    }

    template<typename C1, typename C2>
    MatcherBase<C1, C2>::MatcherBase( const edm::ParameterSet & cfg ) :
      src_( cfg.template getParameter<edm::InputTag>( "src" ) ),
      matched_( cfg.template getParameter<edm::InputTag>( "matched" ) ), 
      distMin_( cfg.template getParameter<double>( "distMin" ) ) {
      produces<MatchMap>();
    }

    template<typename C1, typename C2>
    MatcherBase<C1, C2>::~MatcherBase() { }

    template<typename C1, typename C2>    
    void MatcherBase<C1, C2>::produce( edm::Event& evt, const edm::EventSetup& ) {
      using namespace edm;
      using namespace std;
      Handle<C2> matched;  
      evt.getByLabel( matched_, matched );
      Handle<C1> cands;  
      evt.getByLabel( src_, cands );
      auto_ptr<MatchMap> matchMap( new MatchMap( typename MatchMap::ref_type( RefProd<C1>( cands ), 
									      RefProd<C2>( matched ) ) ) );
      for( size_t c = 0; c != cands->size(); ++ c ) {
	const T1 & cand = (*cands)[ c ];
	vector<helper::MatchPair> v;
	for( size_t m = 0; m != matched->size(); ++ m ) {
	  const T2 & match = ( * matched )[ m ];
	  if ( select( cand, match ) ) {
	    double dist = matchDistance( cand, match );
	    if ( dist < distMin_ ) v.push_back( make_pair( m, dist ) );
	  }
	}
	if ( v.size() > 0 ) {
	  size_t mMin = min_element( v.begin(), v.end(), helper::SortBySecond() )->first;
	  matchMap->insert( Ref<C1>( cands, c ), Ref<C2>( matched, mMin ) );
	}
      }
      
      evt.put( matchMap );
    }    
    
  }
}

#endif
