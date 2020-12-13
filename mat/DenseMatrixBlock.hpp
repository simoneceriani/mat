#pragma once
#include "DenseMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"

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


}