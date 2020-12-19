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

  template< class T, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::DiagonalMatrixBlockIterable() : DiagonalMatrixBlock() {

  }

  // with block
  template< class T, int BR, int BC, int NBR, int NBC >
  template<int Ordering>
  DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::DiagonalMatrixBlockIterable(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp) :
    DiagonalMatrixBlock(blockDesc) 
  {
    // check sparse pattern is diagonal
    for (int o = 0; o < sp.outerSize(); o++) {
      const auto & ins = sp.inner(o);
      assert(ins.size() == 1);
      assert(*(ins.begin()) == o);
    }
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::~DiagonalMatrixBlockIterable() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  template<int Ordering>
  void DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp)  {
    DiagonalMatrixBlock::resize(blockDesc);
    // check sparse pattern is diagonal
    for (int o = 0; o < sp.outerSize(); o++) {
      const auto& ins = sp.inner(o);
      assert(ins.size() == 1);
      assert(*(ins.begin()) == o);
    }
  }


  //----------------------------------------------------------------------------------------------------------------

  template< class T, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT& sm, int id) :
    _curId(id), _sm(&sm), _lastId(id + 1) {

  }
  template< class T, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  DiagonalMatrixBlockIterable<T, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator() {

  }


}