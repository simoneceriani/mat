#include "SparsityPattern.h"

namespace mat {

  template<int Ordering>
  SparsityPattern<Ordering>::SparsityPattern(int nr, int nc) 
    : _count(0)
  {
    if (Ordering == mat::ColMajor) {
      _sp.resize(nc);
      _innerSize = nr;
    }
    else if (Ordering == mat::RowMajor) {
      _sp.resize(nr);
      _innerSize = nc;
    }
    else {
      __MAT_ASSERT_FALSE();
    }
  }

  template<int Ordering>
  SparsityPattern<Ordering>::~SparsityPattern() {

  }

  template<int Ordering>
  Eigen::SparseMatrix<int, Ordering> SparsityPattern<Ordering>::toSparseMatrix() const {

    Eigen::SparseMatrix<int, Ordering> mat;
    if (Ordering == mat::ColMajor) {
      mat.resize(_innerSize, _sp.size());
    }
    else if (Ordering == mat::RowMajor) {
      mat.resize(_sp.size(), _innerSize);
    }
    else {
      __MAT_ASSERT_FALSE();
    }

    mat.reserve(_count);
    for (int o = 0; o < this->outerSize(); o++) {
      mat.startVec(o);
      for (int i : this->inner(o)) {
        mat.insertBackByOuterInner(o, i) = 1; // init with 0
      }
    }
    mat.finalize();
    return mat;
  }
  
  // static
  template<int Ordering>
  typename SparsityPattern<Ordering>::CSPtr SparsityPattern<Ordering>::makeDiagonal(int nr, int nc) {
    auto ret = std::make_shared< SparsityPattern<Ordering>>(nr,nc);
    ret->setDiagonal();
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