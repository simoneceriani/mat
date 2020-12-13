#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"

namespace mat {

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DenseMatrixBlock final : public MatrixBlockBase<BR, BC, NBR, NBC>  {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockDense, T, BR, BC, NBR, NBC>;

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


  public:
    DenseMatrixBlock();

    // with block
    DenseMatrixBlock(const BlockDescriptor& blockDesc);

    virtual ~DenseMatrixBlock();

    void resize(const BlockDescriptor& blockDesc) override final;


    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    inline BlockType block(int r, int c) {
      return _mat.block<RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>(
        this->rowBlockStart(r), this->colBlockStart(r),
        this->rowBlockSize(r), this->colBlockSize(r)
        );
    }

    inline ConstBlockType block(int r, int c) const {
      return _mat.block<RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>(
        this->rowBlockStart(r), this->colBlockStart(r),
        this->rowBlockSize(r), this->colBlockSize(r)
        );
    }

    void setZero();

  };


}