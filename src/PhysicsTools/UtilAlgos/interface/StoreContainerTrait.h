#ifndef UtilAlgos_StoreContainerTrait_h
#define UtilAlgos_StoreContainerTrait_h
/* \class helper::StoreContainerTrait
 *
 * \author Luca Lista, INFN
 */
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"

namespace helper {
  template<typename OutputCollection>
    struct StoreContainerTrait {
      typedef std::vector<const typename OutputCollection::value_type *> type;
  };

  template<typename C>
  struct StoreContainerTrait<edm::RefVector<C> > {
    typedef edm::RefVector<C> type;
  };

  template<typename T>
  struct StoreContainerTrait<std::vector<edm::RefToBase<T> > > {
    typedef std::vector<edm::RefToBase<T> > type;
  };

  template<typename R, typename C>
   struct StoreContainerTrait<edm::AssociationVector<R, C> > {
     typedef typename StoreContainerTrait<typename R::product_type>::type type;
  };
}

#endif
