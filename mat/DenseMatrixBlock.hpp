#pragma once
#include "DenseMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"

namespace mat {

  template< class T, int BR, int BC, int NBR , int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR , int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes)
    : _blockDescriptor(rowBlocksSizes, colBlocksSizes)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR , int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol)
    : _blockDescriptor(rowBlocksSizes, nBlocksRow, colBlocksSizes, nBlocksCol)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::~DenseMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol) {
    _blockDescriptor.resize(rowBlocksSizes, nBlocksRow, colBlocksSizes, nBlocksCol);
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::setZero() {
    _mat.setZero();
  }


}