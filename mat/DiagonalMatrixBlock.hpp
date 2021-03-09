#pragma once
#include "DiagonalMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"

namespace mat {

  template< class T, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::DiagonalMatrixBlock() {
    this->resizeImpl();
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::DiagonalMatrixBlock(const BlockDescriptor& blockDesc)
    : MatrixBlockBase<BR, BC, NBR, NBC>(blockDesc)
  {
    this->resizeImpl();
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::resizeImpl() {
    _mat.resize(std::min(this->numBlocksRow(), this->numBlocksCol()));
    if (RowTraits::blockSizeAtCompileTime == mat::Dynamic || ColTraits::blockSizeAtCompileTime == mat::Dynamic) {
      for (int i = 0; i < _mat.size(); i++) {
        _mat[i].resize(this->rowBlockSize(i), this->colBlockSize(i));
      }
    }
  }


  template< class T, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::~DiagonalMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc) {
    MatrixBlockBase<BR, BC, NBR, NBC>::resize(blockDesc);
    this->resizeImpl();
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DiagonalMatrixBlock<T, BR, BC, NBR, NBC>::setZero() {
    for (int i = 0; i < _mat.size(); i++) {
      _mat[i].setZero();
    }
  }
  //----------------------------------------------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::DiagonalMatrixBlockIterable() : DiagonalMatrixBlockT() {
    auto sp = std::make_shared<SparsityPattern<Ordering>>(this->numBlocksRow(), this->numBlocksCol());
    sp->setDiagonal();
    this->_sparsityPattern = sp;
  }

  // with block
  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::DiagonalMatrixBlockIterable(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp) :
    DiagonalMatrixBlockT(blockDesc),
    _sparsityPattern(sp)
  {
    // check sparse pattern is diagonal
    for (int o = 0; o < std::min(sp->outerSize(), sp->innerSize()); o++) {
      const auto& ins = sp->inner(o);
      assert(ins.size() == 1);
      assert(*(ins.begin()) == o);
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::~DiagonalMatrixBlockIterable() {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp) {
    this->_sparsityPattern = sp;
    DiagonalMatrixBlockT::resize(blockDesc);
    // check sparse pattern is diagonal
    for (int o = 0; o < std::min(sp->outerSize(), sp->innerSize()); o++) {
      const auto& ins = sp->inner(o);
      assert(ins.size() == 1);
      assert(*(ins.begin()) == o);
    }
  }


  //----------------------------------------------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT& sm, int id, int lastId) :
    _curId(id), _sm(&sm), _lastId(lastId) {

  }
  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DiagonalMatrixBlockIterable<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator() {

  }


}