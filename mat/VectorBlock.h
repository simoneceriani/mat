#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockTypeTraits.h"
#include "DimensionDescriptor.h"

namespace mat {

  template<class T, int BR, int NBR = mat::Dynamic >
  class VectorBlock {

  public:

    using Traits = VectorBlockTraits<T, BR, NBR>;

    using StorageType = typename Traits::StorageType;
    using SubVectorType = typename Traits::SubVectorType;
    using SegmentType = typename Traits::SegmentType;
    using ConstSegmentType = typename Traits::ConstSegmentType;

    using SegmentDescriptor = DimensionDescriptor<BR, NBR>;

    using SegmentTraits = typename DimensionDescriptor<BR, NBR>::Traits;

  private:
    StorageType _mat;
    std::shared_ptr<const SegmentDescriptor> _segmentDesc;

  public:
    VectorBlock();

    // with block
    VectorBlock(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes);
    VectorBlock(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow);
    VectorBlock(const std::shared_ptr<const SegmentDescriptor>& rowDescription);

    virtual ~VectorBlock();

    void resize(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow);
    void resize(const std::shared_ptr<const SegmentDescriptor>& rowDescription);


    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    inline SegmentType segment(int r) {
      return _mat.template segment<SegmentTraits::blockSizeAtCompileTime>(
        this->segmentStart(r), this->segmentSize(r)
        );
    }

    inline ConstSegmentType segment(int r) const {
      return _mat.template segment<SegmentTraits::blockSizeAtCompileTime>(
        this->segmentStart(r), this->segmentSize(r)
        );
    }

    void setZero();

    inline int numSegments() const {
      return _segmentDesc->numBlocks();
    }

    inline int segmentUniqueBlockSize() const {
      return _segmentDesc->uniqueBlockSize();
    }

    inline int segmentSize(int i) const {
      return _segmentDesc->blockSize(i);
    }

    inline int segmentStart(int i) const {
      return _segmentDesc->blockStart(i);
    }

    inline int numElements() const {
      return _segmentDesc->numElements();
    }

    const SegmentDescriptor& segmentDescription() const { return *_segmentDesc; }
    const std::shared_ptr<const SegmentDescriptor>& segmentDescriptionCSPtr() const { return _segmentDesc; }

  };


}