#pragma once
#include "SparseCoeffMatrixBlock.h"

#include "SparsityPatternBlockDescriptor.hpp"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"

#include "Utils.hpp"

namespace mat {

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseCoeffMatrixBlock()
  {
    this->_sparsityPattern = std::make_shared<SparsityPattern<Ordering>>(this->numBlocksRow(), this->numBlocksCol());
    _sparseCoeffMap.reset(new SparsityPatternBlockDescriptor<Ordering>(*_sparsityPattern, this->blockDescriptor()));
    this->resizeImpl(*_sparsityPattern);
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseCoeffMatrixBlock(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr & sp)
    : MatrixBlockBase<BR, BC, NBR, NBC>(blockDesc),
    _sparsityPattern(sp)
  {
    if (blockDesc.numBlocksRow() != 0 && blockDesc.numBlocksCol() != 0) {
      _sparseCoeffMap.reset(new SparsityPatternBlockDescriptor<Ordering>(*sp, this->blockDescriptor()));
      this->resizeImpl(*sp);
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::~SparseCoeffMatrixBlock() {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template <class OuterDesc, class InnerDesc>
  void SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::populateSparseMat(const SparsityPattern<Ordering>& sp, const OuterDesc& outDesc, const InnerDesc& inDesc) {
    _mat.resize(this->numElementRows(), this->numElementCols());
    _mat.reserve(_sparseCoeffMap->nonZeroCoeffs());

    // Comments are relative to colmajor
    for (int o = 0; o < sp.outerSize(); o++) {
      // o is the col index
      const int o_start = outDesc.blockStart(o);

      // iterate on the colums of this block, oi is the number of col in the block
      for (int oi = 0; oi < outDesc.blockSize(o); oi++) {

        _mat.startVec(o_start + oi);

        // iterate on blocks in this col
        for (int in : sp.inner(o)) {
          // in is the row index
          const int in_start = inDesc.blockStart(in);

          // iterate on rows of this block
          for (int ini = 0; ini < inDesc.blockSize(in); ini++) {
            // add element to sparse matrix
            _mat.insertBackByOuterInner(o_start + oi, in_start + ini) = 0; // init with 0
          }

        }

      }
    }
    _mat.finalize();

  }


  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resizeImpl(const SparsityPattern<Ordering>& sp) {

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

    // init sparse mat elements

    if (Ordering == mat::ColMajor) {
      populateSparseMat(sp, this->blockDescriptor().colDescription(), this->blockDescriptor().rowDescription());
    }
    else if (Ordering == mat::RowMajor) {
      populateSparseMat(sp, this->blockDescriptor().rowDescription(), this->blockDescriptor().colDescription());
    }
    else __MAT_ASSERT_FALSE();

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr & sp) {
    MatrixBlockBase<BR, BC, NBR, NBC>::resize(blockDesc);
    _sparsityPattern = sp;
    _sparseCoeffMap.reset(new SparsityPatternBlockDescriptor<Ordering>(*sp, this->blockDescriptor()));
    
    this->resizeImpl(*sp);
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::setZero() {
    _mat.coeffs().setZero();
  }
  
  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  int SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::blockUID(int r, int c) const {
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

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  int SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::searchBlockInner(int r, int c) const {
    if (Ordering == mat::RowMajor) {
      int start = _outerStarts[r];
      int end = _outerStarts[r + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = Utils::binary_search(itS, itE, c);
      if (id == itE) return -1;
      else return id - itS;
    }
    else if (Ordering == mat::ColMajor) {
      int start = _outerStarts[c];
      int end = _outerStarts[c + 1];

      auto itS = _innerIndexes.begin() + start;
      auto itE = _innerIndexes.begin() + end;

      auto id = Utils::binary_search(itS, itE, r);
      if (id == itE) return -1;
      else return id - itS;
    }
    else {
      __MAT_ASSERT_FALSE();
    }
    return -1;
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT& sm, int outer, int curId, int lastId)
    : _outer(outer), _sm(&sm), _curId(curId), _lastId(lastId)
  {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseCoeffMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator()
  {

  }


}