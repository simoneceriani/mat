#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>
#include <memory>

#include "DimensionDescriptor.h"

namespace mat {

  template<int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic>
  class MatrixBlockDescriptor {
  public:

    using DimensionDescriptorRow = DimensionDescriptor<BR, NBR>;
    using DimensionDescriptorCol = DimensionDescriptor<BC, NBC>;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

  private:
    
    std::shared_ptr<const DimensionDescriptorRow> _rowDesc;
    std::shared_ptr<const DimensionDescriptorCol> _colDesc;
    
  public:
    MatrixBlockDescriptor();

    // standard
    MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename ColTraits::BlockSizeTypePar colBlocksSizes);
    MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename ColTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);
    MatrixBlockDescriptor(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription);

    // square matrices
    static MatrixBlockDescriptor squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes);
    static MatrixBlockDescriptor squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes, int nBlocks);
    static MatrixBlockDescriptor squareMatrix(const std::shared_ptr<const DimensionDescriptorRow>& blocksDescription);

    virtual ~MatrixBlockDescriptor();

    inline int numBlocksRow() const {
      return _rowDesc->numBlocks();
    }
    inline int numBlocksCol() const {
      return _colDesc->numBlocks();
    }

    inline int rowUniqueBlockSize() const {
      return _rowDesc->uniqueBlockSize();
    }

    inline int colUniqueBlockSize() const {
      return _colDesc->uniqueBlockSize();
    }

    inline int rowBlockSize(int i) const {
      return _rowDesc->blockSize(i);
    }
    inline int colBlockSize(int j) const {
      return _colDesc->blockSize(j);
    }

    inline int rowBlockStart(int i) const {
      return _rowDesc->blockStart(i);
    }
    inline int colBlockStart(int j) const {
      return _colDesc->blockStart(j);
    }

    inline int numElementRows() const {
      return _rowDesc->numElements();
    }
    inline int numElementCols() const {
      return _colDesc->numElements();
    }

    const DimensionDescriptorRow & rowDescription() const { return *_rowDesc; }
    const DimensionDescriptorCol & colDescription() const { return *_colDesc; }

    const std::shared_ptr<const DimensionDescriptorRow>& rowDescriptionCSPtr() const { return _rowDesc; }
    const std::shared_ptr<const DimensionDescriptorCol>& colDescriptionCSPtr() const { return _colDesc; }

    void resizeSquare(const std::shared_ptr<const DimensionDescriptorRow>& desc);
    void resizeSquare(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow);

    void resize(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription);
    void resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename ColTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);

  };

}
