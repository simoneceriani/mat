#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockDescriptor.h"

namespace mat {

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DenseMatrixBlock {

  public:
    using BlockDescriptor = MatrixBlockDescriptor<BR, BC, NBR, NBC>;
    using DimensionDescriptorRow = typename BlockDescriptor::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename BlockDescriptor::DimensionDescriptorCol;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

    using MatrixType = Eigen::Matrix<T, RowTraits::numElementsAtCompileTime, ColTraits::numElementsAtCompileTime>;

  private:
    BlockDescriptor _blockDescriptor;
    MatrixType _mat;

  public:
    DenseMatrixBlock();

    // standard
    DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes);
    DenseMatrixBlock(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);
    DenseMatrixBlock(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription);

    // with block
    DenseMatrixBlock(const BlockDescriptor& blockDesc);

    // square matrices
    static DenseMatrixBlock squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes);
    static DenseMatrixBlock squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes, int nBlocks);
    static DenseMatrixBlock squareMatrix(const std::shared_ptr<const DimensionDescriptorRow>& blocksDescription);


    virtual ~DenseMatrixBlock();

    inline int numBlocksRow() const {
      return _blockDescriptor.numBlocksRow();
    }
    inline int numBlocksCol() const {
      return _blockDescriptor.numBlocksCol();
    }

    inline int rowBlockSize(int i) const {
      return _blockDescriptor.rowBlockSize(i);
    }
    inline int colBlockSize(int j) const {
      return _blockDescriptor.colBlockSize(j);
    }

    inline int rowBlockStart(int i) const {
      return _blockDescriptor.rowBlockStart(i);
    }
    inline int colBlockStart(int j) const {
      return _blockDescriptor.colBlockStart(j);
    }

    inline int numElementRows() const {
      return _blockDescriptor.numElementRows();
    }
    inline int numElementCols() const {
      return _blockDescriptor.numElementCols();
    }

    const BlockDescriptor& blockDescriptor() const {
      return _blockDescriptor;
    }

    void resizeSquare(const std::shared_ptr<const DimensionDescriptorRow>& desc);
    void resizeSquare(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow);

    void resize(const BlockDescriptor& blockDesc);
    void resize(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription);
    void resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);

    const MatrixType& mat() const {
      return _mat;
    }

    void setZero();

  };


}