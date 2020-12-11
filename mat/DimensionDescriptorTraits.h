#pragma once
#include "Global.h"
#include <vector>

namespace mat {

  template<int B, int NB>
  struct DimensionDescriptorTraits {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = B;
    static constexpr int numElementsAtCompileTime = NB * B;
  };

  template<int B>
  struct DimensionDescriptorTraits<B, mat::Dynamic> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = B;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<int NB>
  struct DimensionDescriptorTraits<mat::Dynamic, NB> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = mat::Dynamic;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<>
  struct DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = mat::Dynamic;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<int NB>
  struct DimensionDescriptorTraits<mat::Variable, NB> {
    using BlockSizeType = std::vector<int>;
    using BlockSizeTypePar = const std::vector<int>&;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = mat::Variable;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<>
  struct DimensionDescriptorTraits<mat::Variable, mat::Dynamic> {
    using BlockSizeType = std::vector<int>;
    using BlockSizeTypePar = const std::vector<int>&;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = mat::Variable;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

}