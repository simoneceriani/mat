#pragma once
#include "DenseMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"
#include "Utils.hpp"

namespace mat {

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock() {
    _mat.resize(this->numElementRows(), this->numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(const BlockDescriptor& blockDesc)
    : MatrixBlockBase<BR, BC, NBR, NBC>(blockDesc)
  {
    _mat.resize(this->numElementRows(), this->numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::~DenseMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc) {
    MatrixBlockBase<BR, BC, NBR, NBC>::resize(blockDesc);
    _mat.resize(this->numElementRows(), this->numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::setZero() {
    _mat.setZero();
  }


  //--------------------------------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::DenseMatrixBlockIterable()
    : DenseMatrixBlock<T, BR, BC, NBR, NBC>()
  {
    this->_sparsityPattern = std::make_shared<SparsityPattern<Ordering>>(this->numBlocksRow(), this->numBlocksCol());
    this->createPattern(*_sparsityPattern);
  }


  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::DenseMatrixBlockIterable(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp)
    : DenseMatrixBlock <T, BR, BC, NBR, NBC>(blockDesc)
  {
    this->_sparsityPattern = sp;
    if (blockDesc.numBlocksRow() != 0 && blockDesc.numBlocksCol() != 0) {
      this->createPattern(*sp);
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::~DenseMatrixBlockIterable() {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp) {
    DenseMatrixBlock<T, BR, BC, NBR, NBC>::resize(blockDesc);
    this->_sparsityPattern = sp;
    this->createPattern(*sp);
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::createPattern(const SparsityPattern<Ordering>& sp) {
    // init block structure
    _innerIndexes.resize(sp.count());
    _uid2outer.resize(sp.count());
    _outerStarts.resize(sp.outerSize() + 1);
    int count = 0;
    for (int o = 0; o < sp.outerSize(); o++) {
      _outerStarts[o] = count;
      for (int in : sp.inner(o)) {
        _innerIndexes[count] = in;
        _uid2outer[count] = o;
        count++;
      }
    }
    _outerStarts[sp.outerSize()] = count;
    assert(count == sp.count());

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  int DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::blockUID(int r, int c) const {
    if (Ordering == mat::RowMajor) {
      int start = _outerStarts[r];
      int end = _outerStarts[r + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = Utils::binary_search(itS, itE, c);
      if (id == itE) return -1;
      else return id - _innerIndexes.begin();
    }
    else if (Ordering == mat::ColMajor) {
      int start = _outerStarts[c];
      int end = _outerStarts[c + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = Utils::binary_search(itS, itE, r);
      if (id == itE) return -1;
      else return id - _innerIndexes.begin();
    }
    else {
      __MAT_ASSERT_FALSE();
    }
    return -1;
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT& sm, int outer, int curId, int lastId)
    : _outer(outer), _sm(&sm), _curId(curId), _lastId(lastId)
  {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DenseMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator()
  {

  }


}