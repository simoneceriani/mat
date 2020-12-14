#pragma once
#include "SparseMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"

namespace mat {

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseMatrixBlock() {
    this->resizeImpl(SparsityPattern<Ordering>(this->numBlocksRow(), this->numBlocksCol()));
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseMatrixBlock(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp)
    : MatrixBlockBase<BR, BC, NBR, NBC>(blockDesc)
  {
    if (blockDesc.numBlocksRow() != 0 && blockDesc.numBlocksCol() != 0) {
      this->resizeImpl(sp);
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::~SparseMatrixBlock() {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resizeImpl(const SparsityPattern<Ordering>& sp) {
    _mat.resize(sp.count());
    _innerIndexes.resize(sp.count());
    _outerStarts.resize(sp.outerSize() + 1);

    int count = 0;
    for (int o = 0; o < sp.outerSize(); o++) {
      _outerStarts[o] = count;
      for (int in : sp.inner(o)) {
        _innerIndexes[count] = in;

        if (RowTraits::blockSizeAtCompileTime == mat::Dynamic || ColTraits::blockSizeAtCompileTime == mat::Dynamic) {
          if (Ordering == mat::RowMajor) {
            _mat[count].resize(this->rowBlockSize(o), this->colBlockSize(in));
          }
          else if (Ordering == mat::ColMajor) {
            _mat[count].resize(this->rowBlockSize(in), this->colBlockSize(o));
          }
          else {
            assert(false && "how do you end up here?!?");
          }
        }

        count++;
      }
    }
    _outerStarts[sp.outerSize()] = count;
    assert(count == sp.count());

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp) {
    MatrixBlockBase<BR, BC, NBR, NBC>::resize(blockDesc);
    this->resizeImpl(sp);
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::setZero() {
    for (int i = 0; i < _mat.size(); i++) {
      _mat[i].setZero();
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template <class ForwardIterator>
  static ForwardIterator SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::binary_search(ForwardIterator first, ForwardIterator last, int val)
  {
    first = std::lower_bound(first, last, val);
    if (first != last) {
      if (val == *first) return first;
      return last;
    }
    return last;
  }


  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  int SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::searchBlockUID(int r, int c) {
    if (Ordering == mat::RowMajor) {
      int start = _outerStarts[r];
      int end = _outerStarts[r + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = binary_search(itS, itE, c);
      if (id == itE) return -1;
      else return id - _innerIndexes.begin();
    }
    else if (Ordering == mat::ColMajor) {
      int start = _outerStarts[c];
      int end = _outerStarts[c + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = binary_search(itS, itE, r);
      if (id == itE) return -1;
      else return id - _innerIndexes.begin();
    }
    else {
      assert(false && "how do you end up here?!?");
    }
    return -1;
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT & sm, int curId, int lastId)
    : _sm(&sm), _curId(curId), _lastId(lastId) 
  {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator()
  {

  }


}