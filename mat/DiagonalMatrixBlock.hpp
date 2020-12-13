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
    _mat.resize(std::min(numBlocksRow(), numBlocksCol()));
    if (RowTraits::blockSizeAtCompileTime == mat::Dynamic || ColTraits::blockSizeAtCompileTime == mat::Dynamic) {
      for (int i = 0; i < _mat.size(); i++) {
        _mat[i].resize(rowBlockSize(i), colBlockSize(i));
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
      _mat.setZero();
    }    
  }

}