#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"

namespace mat {

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DiagonalMatrixBlock final : public MatrixBlockBase<BR, BC, NBR, NBC>  {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockDiagonal, T, BR, BC, NBR, NBC>;

    using StorageType = typename Traits::StorageType;
    using SubMatrixType = typename Traits::SubMatrixType;
    using BlockType = typename Traits::BlockType;
    using ConstBlockType = typename Traits::ConstBlockType;

    using BlockDescriptor = MatrixBlockDescriptor<BR, BC, NBR, NBC>;
    using DimensionDescriptorRow = typename BlockDescriptor::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename BlockDescriptor::DimensionDescriptorCol;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

  private:
    StorageType _mat;

    void resizeImpl();

  public:
    DiagonalMatrixBlock();

    // with block
    DiagonalMatrixBlock(const BlockDescriptor& blockDesc);

    virtual ~DiagonalMatrixBlock();

    void resize(const BlockDescriptor& blockDesc) override final;


    const StorageType & mat() const {
      return _mat;
    }

    StorageType & mat() {
      return _mat;
    }

    inline BlockType block(int r, int c) {
      assert(r == c);
      assert(r >= 0 && r < _mat.size());
      return _mat[r];
    }

    inline ConstBlockType block(int r, int c) const {
      assert(r == c);
      assert(r >= 0 && r < _mat.size());
      return _mat[r];
    }

    void setZero();

  };


}