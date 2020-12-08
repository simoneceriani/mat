#include "DimensionDescriptor.h"
#include "DimensionDescriptor.hpp"

namespace mat {

  DimensionDescriptorBase<mat::Dynamic>::DimensionDescriptorBase() : _nb(0) {

  }

  DimensionDescriptorBase<mat::Dynamic>::DimensionDescriptorBase(int nb) : _nb(nb) {

  }

  DimensionDescriptorBase<mat::Dynamic>::~DimensionDescriptorBase() {

  }

  void DimensionDescriptorBase<mat::Dynamic>::setNumBlocks(int nb) {
    this->_nb = nb;
  }


}