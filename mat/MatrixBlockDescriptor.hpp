#pragma once
#include "MatrixBlockDescriptor.h"

#include "DimensionDescriptor.hpp"

namespace mat {

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor()
    : _rowDesc(std::make_shared<const DimensionDescriptorRow>()), _colDesc(std::make_shared <const DimensionDescriptorCol>())
  {

  }

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes)
    : _rowDesc(std::make_shared<const DimensionDescriptorRow>(rowBlocksSizes)), _colDesc(std::make_shared <const DimensionDescriptorCol>(colBlocksSizes))
  {

  }

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol)
    : _rowDesc(std::make_shared<const DimensionDescriptorRow>(rowBlocksSizes, nBlocksRow)), _colDesc(std::make_shared <const DimensionDescriptorCol>(colBlocksSizes, nBlocksCol))
  {

  }

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription)
    : _rowDesc(rowDescription), _colDesc(colDescription)
  {

  }

  // square matrices
  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC> MatrixBlockDescriptor<BR, BC, NBR, NBC>::squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes) 
  {
    auto desc = std::make_shared<const DimensionDescriptorRow>(blocksSizes);
    return MatrixBlockDescriptor<BR, BC, NBR, NBC>(desc, desc);
  }

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC> MatrixBlockDescriptor<BR, BC, NBR, NBC>::squareMatrix(typename RowTraits::BlockSizeTypePar blocksSizes, int nBlocks)
  {
    auto desc = std::make_shared<const DimensionDescriptorRow>(blocksSizes, nBlocks);
    return MatrixBlockDescriptor<BR, BC, NBR, NBC>(desc, desc);
  }


  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC> MatrixBlockDescriptor<BR, BC, NBR, NBC>::squareMatrix(const std::shared_ptr<const DimensionDescriptorRow>& blocksDescription)
  {
    return MatrixBlockDescriptor<BR, BC, NBR, NBC>(blocksDescription, blocksDescription);
  }


  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::~MatrixBlockDescriptor() {

  }

  template<int BR, int BC, int NBR, int NBC>
  void MatrixBlockDescriptor<BR, BC, NBR, NBC>::resizeSquare(const std::shared_ptr<const DimensionDescriptorRow>& desc) {
    this->resize(desc, desc);
  }
  
  template<int BR, int BC, int NBR, int NBC>
  void MatrixBlockDescriptor<BR, BC, NBR, NBC>::resizeSquare(typename RowTraits::BlockSizeTypePar blocksSizes, int nBlocks) {
    _rowDesc = std::make_shared<const DimensionDescriptorRow>(blocksSizes, nBlocks);
    _colDesc = _rowDesc;
  }

  template<int BR, int BC, int NBR, int NBC>
  void MatrixBlockDescriptor<BR, BC, NBR, NBC>::resize(const std::shared_ptr<const DimensionDescriptorRow>& rowDescription, const std::shared_ptr<const DimensionDescriptorCol>& colDescription) {
    _rowDesc = rowDescription;
    _colDesc = colDescription;
  }

  template<int BR, int BC, int NBR, int NBC>
  void MatrixBlockDescriptor<BR, BC, NBR, NBC>::resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol) {
    _rowDesc = std::make_shared<const DimensionDescriptorRow>(rowBlocksSizes, nBlocksRow);
    _colDesc = std::make_shared<const DimensionDescriptorCol>(colBlocksSizes, nBlocksCol);
  }

}
