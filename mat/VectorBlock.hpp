#pragma once
#include "VectorBlock.h"

#include "DimensionDescriptor.hpp"

namespace mat {

  template<class T, int BR, int NBR>
  VectorBlock<T, BR, NBR>::VectorBlock() :
    _segmentDesc(std::make_shared<const SegmentDescriptor>())
  {
    _mat.resize(this->numElements());
  }

  template<class T, int BR, int NBR>
  VectorBlock<T, BR, NBR>::VectorBlock(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes) :
    _segmentDesc(std::make_shared<const SegmentDescriptor>(rowBlocksSizes))

  {
    _mat.resize(this->numElements());
  }

  template<class T, int BR, int NBR>
  VectorBlock<T, BR, NBR>::VectorBlock(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow) :
    _segmentDesc(std::make_shared<const SegmentDescriptor>(rowBlocksSizes, nBlocksRow))
  {
    _mat.resize(this->numElements());
  }

  template<class T, int BR, int NBR>
  VectorBlock<T, BR, NBR>::VectorBlock(const std::shared_ptr<const SegmentDescriptor>& rowDescription) :
    _segmentDesc(rowDescription)
  {
    _mat.resize(this->numElements());
  }

  //-----------------------------------------------------------

  template<class T, int BR, int NBR>
  VectorBlock<T, BR, NBR>::~VectorBlock()
  {

  }

  //-----------------------------------------------------------

  template<class T, int BR, int NBR>
  void VectorBlock<T, BR, NBR>::resize(typename SegmentTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow) {
    _segmentDesc = std::make_shared<const SegmentDescriptor>(rowBlocksSizes, nBlocksRow);
    _mat.resize(this->numElements());
  }

  template<class T, int BR, int NBR>
  void VectorBlock<T, BR, NBR>::resize(const std::shared_ptr<const SegmentDescriptor>& rowDescription) {
    _segmentDesc = rowDescription;
    _mat.resize(this->numElements());
  }

  template<class T, int BR, int NBR>
  void VectorBlock<T, BR, NBR>::setZero() {
    _mat.setZero();
  }


}