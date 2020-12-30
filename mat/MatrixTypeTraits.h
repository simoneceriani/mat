#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include "DenseMatrixBlock.h"
#include "DiagonalMatrixBlock.h"
#include "SparseCoeffMatrixBlock.h"
#include "SparseMatrixBlock.h"

namespace mat {

  template<int matType, class T, int Ordering, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic>
  class MatrixBlockIterableTypeTraits {

  };

  template<class T, int Ordering, int BR, int BC, int NBR, int NBC>
  class MatrixBlockIterableTypeTraits<mat::BlockDense, T, Ordering, BR, BC, NBR, NBC> {
  public:
    using MatrixType = DenseMatrixBlockIterable < T, Ordering, BR, BC, NBR, NBC>;
  };

  template<class T, int Ordering, int BR, int BC, int NBR, int NBC>
  class MatrixBlockIterableTypeTraits<mat::BlockDiagonal, T, Ordering, BR, BC, NBR, NBC> {
  public:
    using MatrixType = DiagonalMatrixBlockIterable< T, Ordering, BR, BC, NBR, NBC>;
  };

  template<class T, int Ordering, int BR, int BC, int NBR, int NBC>
  class MatrixBlockIterableTypeTraits<mat::BlockSparse, T, Ordering, BR, BC, NBR, NBC> {
  public:
    using MatrixType = SparseMatrixBlock< T, Ordering, BR, BC, NBR, NBC>;
  };

  template<class T, int Ordering, int BR, int BC, int NBR, int NBC>
  class MatrixBlockIterableTypeTraits<mat::BlockCoeffSparse, T, Ordering, BR, BC, NBR, NBC> {
  public:
    using MatrixType = SparseCoeffMatrixBlock< T, Ordering, BR, BC, NBR, NBC>;
  };

}