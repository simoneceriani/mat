#pragma once
#include "DimensionDescriptorTraits.h"

namespace mat {

  template<int B, int NB>
  constexpr int DimensionDescriptorTraits<B,NB>::numBlocksAtCompileTime;
  template<int B, int NB>
  constexpr int DimensionDescriptorTraits<B,NB>::blockSizeAtCompileTime;
  template<int B, int NB>
  constexpr int DimensionDescriptorTraits<B,NB>::numElementsAtCompileTime;

  template<int B>
  constexpr int DimensionDescriptorTraits<B,mat::Dynamic>::numBlocksAtCompileTime;
  template<int B>
  constexpr int DimensionDescriptorTraits<B,mat::Dynamic>::blockSizeAtCompileTime;
  template<int B>
  constexpr int DimensionDescriptorTraits<B,mat::Dynamic>::numElementsAtCompileTime;

  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Dynamic,NB>::numBlocksAtCompileTime;
  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Dynamic,NB>::blockSizeAtCompileTime;
  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Dynamic,NB>::numElementsAtCompileTime;

  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Variable,NB>::numBlocksAtCompileTime;
  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Variable,NB>::blockSizeAtCompileTime;
  template<int NB>
  constexpr int DimensionDescriptorTraits<mat::Variable,NB>::numElementsAtCompileTime;

}