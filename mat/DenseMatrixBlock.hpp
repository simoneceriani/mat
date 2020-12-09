#pragma once
#include "DenseMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"

namespace mat {

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes)
    : _blockDescriptor(rowBlocksSizes, colBlocksSizes)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol)
    : _blockDescriptor(rowBlocksSizes, nBlocksRow, colBlocksSizes, nBlocksCol)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription)
    : _blockDescriptor(rowDescription, colDescription)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::DenseMatrixBlock(const BlockDescriptor& blockDesc)
    : _blockDescriptor(blockDesc)
  {
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC> DenseMatrixBlock<T, BR, BC, NBR, NBC>::squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes) {
    return DenseMatrixBlock<T, BR, BC, NBR, NBC>(BlockDescriptor::squareMatrix(blocksSizes));
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC> DenseMatrixBlock<T, BR, BC, NBR, NBC>::squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes, int nBlocks) {
    return DenseMatrixBlock<T, BR, BC, NBR, NBC>(BlockDescriptor::squareMatrix(blocksSizes, nBlocks));
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC> DenseMatrixBlock<T, BR, BC, NBR, NBC>::squareMatrix(const std::shared_ptr<const DimensionDescriptorRow>& blocksDescription) {
    return DenseMatrixBlock<T, BR, BC, NBR, NBC>(blocksDescription, blocksDescription);
  }



  template< class T, int BR, int BC, int NBR, int NBC >
  DenseMatrixBlock<T, BR, BC, NBR, NBC>::~DenseMatrixBlock() {

  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resizeSquare(const std::shared_ptr<const DimensionDescriptorRow>& desc) {
    _blockDescriptor.resizeSquare(desc);
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resizeSquare(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow) {
    _blockDescriptor.resizeSquare(rowBlocksSizes, nBlocksRow);
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }


  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc) {
    _blockDescriptor = blockDesc;
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
  }

  template< class T, int BR, int BC, int NBR, int NBC >
  void DenseMatrixBlock<T, BR, BC, NBR, NBC>::resize(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription) {
    _blockDescriptor.resize(rowDescription, colDescription);
    _mat.resize(_blockDescriptor.numElementRows(), _blockDescriptor.numElementCols());
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