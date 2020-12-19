#include "SparsityPattern.h"

namespace mat {

  template<int Ordering>
  SparsityPattern<Ordering>::SparsityPattern(int nr, int nc) 
    : _count(0)
  {
    if (Ordering == mat::ColMajor) {
      _sp.resize(nc);
    }
    else if (Ordering == mat::RowMajor) {
      _sp.resize(nr);
    }
    else {
      ASSERT_FALSE();
    }
  }

  template<int Ordering>
  SparsityPattern<Ordering>::~SparsityPattern() {

  }
  
  // static
  template<int Ordering>
  SparsityPattern<Ordering> SparsityPattern<Ordering>::makeDiagonal(int nr, int nc) {
    int n = std::min(nr, nc);
    SparsityPattern<Ordering> ret(n,n);
    for (int i = 0; i < ret.outerSize(); i++) {
      ret.add(i, i);
    }
    return ret;
  }

  template<int Ordering>
  void SparsityPattern<Ordering>::clear() {
    for (int k = 0; k < _sp.size(); k++) {
      _sp[k].clear();
    }
    _count = 0;
  }

  // explicit instantiation
  template class SparsityPattern<mat::ColMajor>;
  template class SparsityPattern<mat::RowMajor>;

}