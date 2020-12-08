#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "DimensionDescriptor.h"

namespace mat {

  template<int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic>
  class MatrixBlockDescriptor {
  public:

    using DimensionDescriptorRow = DimensionDescriptor<BR, NBR>;
    using DimensionDescriptorCol = DimensionDescriptor<BC, NBC>;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorRow::Traits;

  private:
    
    DimensionDescriptorRow _rowDesc;
    DimensionDescriptorCol _colDesc;
    
  public:
    MatrixBlockDescriptor();
    MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes);
    MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);
    virtual ~MatrixBlockDescriptor();

    inline int numBlocksRow() const {
      return _rowDesc.numBlocks();
    }
    inline int numBlocksCol() const {
      return _colDesc.numBlocks();
    }

    inline int rowBlockSize(int i) const {
      return _rowDesc.blockSize(i);
    }
    inline int colBlockSize(int j) const {
      return _colDesc.blockSize(j);
    }

    inline int rowBlockStart(int i) const {
      return _rowDesc.blockStart(i);
    }
    inline int colBlockStart(int j) const {
      return _colDesc.blockStart(j);
    }

    inline int numElementRows() const {
      return _rowDesc.numElements();
    }
    inline int numElementCols() const {
      return _colDesc.numElements();
    }

    void resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol);

  };

}
