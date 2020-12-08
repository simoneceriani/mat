#pragma once
#include "MatrixBlockDescriptor.h"

#include "DimensionDescriptor.hpp"

namespace mat {

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor() 
  {

  }

  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, typename RowTraits::BlockSizeTypePar colBlocksSizes) 
    : _rowDesc(rowBlocksSizes), _colDesc(colBlocksSizes)
  {

  }
  
  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::MatrixBlockDescriptor(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol) 
    : _rowDesc(rowBlocksSizes, nBlocksRow), _colDesc(colBlocksSizes, nBlocksCol)
  {

  }
  
  template<int BR, int BC, int NBR, int NBC>
  MatrixBlockDescriptor<BR, BC, NBR, NBC>::~MatrixBlockDescriptor() {

  }

  template<int BR, int BC, int NBR, int NBC>
  void MatrixBlockDescriptor<BR, BC, NBR, NBC>::resize(typename RowTraits::BlockSizeTypePar rowBlocksSizes, int nBlocksRow, typename RowTraits::BlockSizeTypePar colBlocksSizes, int nBlocksCol) {
    _rowDesc.resize(rowBlocksSizes, nBlocksRow);
    _colDesc.resize(colBlocksSizes, nBlocksCol);
  }

}
