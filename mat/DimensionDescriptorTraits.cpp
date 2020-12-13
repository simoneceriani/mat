#include "DimensionDescriptorTraits.h"
#include "DimensionDescriptorTraits.hpp"

namespace mat {

  constexpr int DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic>::numBlocksAtCompileTime;
  constexpr int DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic>::blockType;
  constexpr int DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic>::blockSizeAtCompileTime;
  constexpr int DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic>::numElementsAtCompileTime;

  constexpr int DimensionDescriptorTraits<mat::Variable, mat::Dynamic>::numBlocksAtCompileTime;
  constexpr int DimensionDescriptorTraits<mat::Variable, mat::Dynamic>::blockType;
  constexpr int DimensionDescriptorTraits<mat::Variable, mat::Dynamic>::blockSizeAtCompileTime;
  constexpr int DimensionDescriptorTraits<mat::Variable, mat::Dynamic>::numElementsAtCompileTime;
    
}
