#include "SparsityPattern.h"

namespace mat {

  template<int Ordering>
  SparsityPattern<Ordering>::SparsityPattern(int nr, int nc) 
    : _count(-1)
  {
    if (Ordering == mat::ColMajor) {
      _sp.resize(nc);
    }
    else if (Ordering == mat::RowMajor) {
      _sp.resize(nr);
    }
    else {
      assert(false && "how do you end up here?!?");
    }
  }

  template<int Ordering>
  SparsityPattern<Ordering>::~SparsityPattern() {

  }

  template<int Ordering>
  void SparsityPattern<Ordering>::clear() {
    for (int k = 0; k < _sp.size(); k++) {
      _sp[k].clear();
    }
    _count = -1;
  }

  // explicit instantiation
  template class SparsityPattern<mat::ColMajor>;
  template class SparsityPattern<mat::RowMajor>;

}