#pragma once
#include "MatrixBlockBase.h"

#include "MatrixBlockDescriptor.hpp"

namespace mat {

  template< int BR, int BC, int NBR, int NBC >
  MatrixBlockBase<BR, BC, NBR, NBC>::MatrixBlockBase() {

  }

  template< int BR, int BC, int NBR, int NBC >
  MatrixBlockBase<BR, BC, NBR, NBC>::MatrixBlockBase(const BlockDescriptor& blockDesc)
    : _blockDescriptor(blockDesc)
  {
  }

  template< int BR, int BC, int NBR, int NBC >
  MatrixBlockBase<BR, BC, NBR, NBC>::~MatrixBlockBase() {

  }

  template< int BR, int BC, int NBR, int NBC >
  void MatrixBlockBase<BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc) {
    _blockDescriptor = blockDesc;
  }



}