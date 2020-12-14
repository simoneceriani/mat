#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

namespace mat {

  template<int matType, class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic>
  class MatrixBlockTypeTraits {

  };

  template<class T, int BR, int BC, int NBR, int NBC>
  class MatrixBlockTypeTraits<mat::BlockDense, T, BR, BC, NBR, NBC> {
  public:
    using RowTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorRow::Traits;
    using ColTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorCol::Traits;


    using StorageType = Eigen::Matrix<T, RowTraits::numElementsAtCompileTime, ColTraits::numElementsAtCompileTime>;
    using SubMatrixType = Eigen::Matrix < T, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>;
    using BlockType = typename Eigen::Block<StorageType, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>;
    using ConstBlockType = typename Eigen::Block<const StorageType, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>;

  };

  template<class T, int BR, int BC, int NBR, int NBC>
  class MatrixBlockTypeTraits<mat::BlockDiagonal, T, BR, BC, NBR, NBC> {
  public:
    using RowTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorRow::Traits;
    using ColTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorCol::Traits;


    using StorageType = std::vector<Eigen::Matrix<T, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>>;
    using SubMatrixType = Eigen::Matrix < T, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>;
    using BlockType = SubMatrixType &;
    using ConstBlockType = const SubMatrixType &;

  };

  template<class T, int BR, int BC, int NBR, int NBC>
  class MatrixBlockTypeTraits<mat::BlockSparse, T, BR, BC, NBR, NBC> {
  public:
    using RowTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorRow::Traits;
    using ColTraits = typename MatrixBlockDescriptor<BR, BC, NBR, NBC>::DimensionDescriptorCol::Traits;


    using StorageType = std::vector<Eigen::Matrix<T, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>>;
    using SubMatrixType = Eigen::Matrix < T, RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>;
    using BlockType = SubMatrixType&;
    using ConstBlockType = const SubMatrixType&;

  };

}