#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/SparsityPattern.h"
#include "mat/SparseMatrixBlock.hpp"

// let use W matrix in BA case
static constexpr int camSize = 6;
static constexpr int pointSize = 3;
static constexpr int numCams = 10;
static constexpr int numPoints = 13;

TEMPLATE_TEST_CASE_SIG("SparseMatrixBlock", "[SparseMatrixBlock]", ((int Ordering, int BR, int BC, int NBR, int NBC), Ordering, BR, BC, NBR, NBC),
  (mat::RowMajor, camSize, pointSize, numCams, numPoints),
  (mat::RowMajor, camSize, pointSize, numCams, mat::Dynamic),
  (mat::RowMajor, camSize, pointSize, mat::Dynamic, numPoints),
  (mat::RowMajor, camSize, pointSize, mat::Dynamic, mat::Dynamic),
  (mat::RowMajor, camSize, mat::Dynamic, numCams, numPoints),
  (mat::RowMajor, camSize, mat::Dynamic, numCams, mat::Dynamic),
  (mat::RowMajor, camSize, mat::Dynamic, mat::Dynamic, numPoints),
  (mat::RowMajor, camSize, mat::Dynamic, mat::Dynamic, mat::Dynamic),
  (mat::RowMajor, mat::Dynamic, pointSize, numCams, numPoints),
  (mat::RowMajor, mat::Dynamic, pointSize, numCams, mat::Dynamic),
  (mat::RowMajor, mat::Dynamic, pointSize, mat::Dynamic, numPoints),
  (mat::RowMajor, mat::Dynamic, pointSize, mat::Dynamic, mat::Dynamic),
  (mat::RowMajor, mat::Dynamic, mat::Dynamic, numCams, numPoints),
  (mat::RowMajor, mat::Dynamic, mat::Dynamic, numCams, mat::Dynamic),
  (mat::RowMajor, mat::Dynamic, mat::Dynamic, mat::Dynamic, numPoints),
  (mat::RowMajor, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::Dynamic),
  (mat::ColMajor, camSize, pointSize, numCams, numPoints),
  (mat::ColMajor, camSize, pointSize, numCams, mat::Dynamic),
  (mat::ColMajor, camSize, pointSize, mat::Dynamic, numPoints),
  (mat::ColMajor, camSize, pointSize, mat::Dynamic, mat::Dynamic),
  (mat::ColMajor, camSize, mat::Dynamic, numCams, numPoints),
  (mat::ColMajor, camSize, mat::Dynamic, numCams, mat::Dynamic),
  (mat::ColMajor, camSize, mat::Dynamic, mat::Dynamic, numPoints),
  (mat::ColMajor, camSize, mat::Dynamic, mat::Dynamic, mat::Dynamic),
  (mat::ColMajor, mat::Dynamic, pointSize, numCams, numPoints),
  (mat::ColMajor, mat::Dynamic, pointSize, numCams, mat::Dynamic),
  (mat::ColMajor, mat::Dynamic, pointSize, mat::Dynamic, numPoints),
  (mat::ColMajor, mat::Dynamic, pointSize, mat::Dynamic, mat::Dynamic),
  (mat::ColMajor, mat::Dynamic, mat::Dynamic, numCams, numPoints),
  (mat::ColMajor, mat::Dynamic, mat::Dynamic, numCams, mat::Dynamic),
  (mat::ColMajor, mat::Dynamic, mat::Dynamic, mat::Dynamic, numPoints),
  (mat::ColMajor, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::Dynamic)
)
{
  mat::SparsityPattern<Ordering> sp(numCams, numPoints);
  int count = 0;
  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < numCams; i++) {
      for (int j = 0; j < numPoints; j++) {
        if (((i + j) % 2) == 0) {
          if (!sp.has(i, j)) count++;
          sp.add(i, j);
          REQUIRE(sp.has(i, j));
        }
        else {
          REQUIRE(!sp.has(i, j));
        }
      }
    }
  }

  REQUIRE(sp.count() == count);

  using MatT = mat::SparseMatrixBlock<double, Ordering, BR, BC, NBR, NBC>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));

    REQUIRE(mat->mat().size() == 0);

    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), sp);
    REQUIRE(mat->mat().size() == sp.count());

  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, pointSize), sp));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), sp);
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), sp));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), sp);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

}


TEMPLATE_TEST_CASE_SIG("SparseMatrixBlock-RowMajor", "[SparseMatrixBlock]", ((int BR, int BC, int NBR, int NBC), BR, BC, NBR, NBC),
  (2, 3, 3, 4),
  (2, 3, 3, Eigen::Dynamic),
  (2, 3, Eigen::Dynamic, 4),
  (2, 3, Eigen::Dynamic, Eigen::Dynamic),
  (2, Eigen::Dynamic, 3, 4),
  (2, Eigen::Dynamic, 3, Eigen::Dynamic),
  (2, Eigen::Dynamic, Eigen::Dynamic, 4),
  (2, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic),
  (Eigen::Dynamic, 3, 3, 4),
  (Eigen::Dynamic, 3, 3, Eigen::Dynamic),
  (Eigen::Dynamic, 3, Eigen::Dynamic, 4),
  (Eigen::Dynamic, 3, Eigen::Dynamic, Eigen::Dynamic),
  (Eigen::Dynamic, Eigen::Dynamic, 3, 4),
  (Eigen::Dynamic, Eigen::Dynamic, 3, Eigen::Dynamic),
  (Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, 4),
  (Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic)
) {
  // +-----------------------+
  // | 0.0 | --- | --- | 0.3 |
  // | 1.0 | --- | 1.2 | 1.3 |
  // | 2.0 | --- | 2.2 | --- |
  // +-----------------------+
  mat::SparsityPattern<mat::RowMajor> sp(3, 4);
  sp.add(0, 0); sp.add(0, 3);
  sp.add(1, 0); sp.add(1, 2); sp.add(1, 3);
  sp.add(2, 0); sp.add(2, 2);

  using MatT = mat::SparseMatrixBlock<double, mat::RowMajor, BR, BC, NBR, NBC>;
  MatT mat(typename MatT::BlockDescriptor(2, 3, 3, 4), sp);

  REQUIRE(mat.nonZeroBlocks() == sp.count());

  REQUIRE(mat.searchBlockUID(0, 0) == 0);
  REQUIRE(mat.searchBlockUID(0, 1) == -1);
  REQUIRE(mat.searchBlockUID(0, 2) == -1);
  REQUIRE(mat.searchBlockUID(0, 3) == 1);

  REQUIRE(mat.searchBlockUID(1, 0) == 2);
  REQUIRE(mat.searchBlockUID(1, 1) == -1);
  REQUIRE(mat.searchBlockUID(1, 2) == 3);
  REQUIRE(mat.searchBlockUID(1, 3) == 4);

  REQUIRE(mat.searchBlockUID(2, 0) == 5);
  REQUIRE(mat.searchBlockUID(2, 1) == -1);
  REQUIRE(mat.searchBlockUID(2, 2) == 6);
  REQUIRE(mat.searchBlockUID(2, 3) == -1);

  // inefficient!
  for (int r = 0; r < mat.numBlocksRow(); r++) {
    for (int c = 0; c < mat.numBlocksCol(); c++) {
      if (mat.hasBlock(r, c)) {
        mat.block(r, c).setConstant(r + c / 1.0);
      }
    }
  }

  for (int r = 0; r < mat.numBlocksRow(); r++) {
    for (int c = 0; c < mat.numBlocksCol(); c++) {
      if (mat.hasBlock(r, c)) {
        REQUIRE((mat.block(r, c).array() == r + c / 1.0).all());
      }
    }
  }

  {
    auto it = mat.rowBegin(0);
    REQUIRE(it() == 0);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 1);
    REQUIRE(it.col() == 3);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = mat.rowBegin(1);
    REQUIRE(it() == 2);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 3);
    REQUIRE(it.col() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 4);
    REQUIRE(it.col() == 3);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = mat.rowBegin(2);
    REQUIRE(it() == 5);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 6);
    REQUIRE(it.col() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());
  }

  {
    const auto& matC = mat;
    auto it = matC.rowBegin(0);
    REQUIRE(it() == 0);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 1);
    REQUIRE(it.col() == 3);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = matC.rowBegin(1);
    REQUIRE(it() == 2);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 3);
    REQUIRE(it.col() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 4);
    REQUIRE(it.col() == 3);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = matC.rowBegin(2);
    REQUIRE(it() == 5);
    REQUIRE(it.col() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 6);
    REQUIRE(it.col() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());
  }

}

//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------

TEMPLATE_TEST_CASE_SIG("SparseMatrixBlock-ColMajor", "[SparseMatrixBlock]", ((int BR, int BC, int NBR, int NBC), BR, BC, NBR, NBC),
  (2, 3, 3, 4)/*,
  (2, 3, 3, Eigen::Dynamic),
  (2, 3, Eigen::Dynamic, 4),
  (2, 3, Eigen::Dynamic, Eigen::Dynamic),
  (2, Eigen::Dynamic, 3, 4),
  (2, Eigen::Dynamic, 3, Eigen::Dynamic),
  (2, Eigen::Dynamic, Eigen::Dynamic, 4),
  (2, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic),
  (Eigen::Dynamic, 3, 3, 4),
  (Eigen::Dynamic, 3, 3, Eigen::Dynamic),
  (Eigen::Dynamic, 3, Eigen::Dynamic, 4),
  (Eigen::Dynamic, 3, Eigen::Dynamic, Eigen::Dynamic),
  (Eigen::Dynamic, Eigen::Dynamic, 3, 4),
  (Eigen::Dynamic, Eigen::Dynamic, 3, Eigen::Dynamic),
  (Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, 4),
  (Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic)*/
) {
  // +-----------------------+
  // | 0.0 | --- | --- | 0.3 |
  // | 1.0 | --- | 1.2 | 1.3 |
  // | 2.0 | --- | 2.2 | --- |
  // +-----------------------+
  mat::SparsityPattern<mat::ColMajor> sp(3, 4);
  sp.add(0, 0); sp.add(0, 3);
  sp.add(1, 0); sp.add(1, 2); sp.add(1, 3);
  sp.add(2, 0); sp.add(2, 2);

  using MatT = mat::SparseMatrixBlock<double, mat::ColMajor, BR, BC, NBR, NBC>;
  MatT mat(typename MatT::BlockDescriptor(2, 3, 3, 4), sp);

  REQUIRE(mat.nonZeroBlocks() == sp.count());

  REQUIRE(mat.searchBlockUID(0, 0) == 0);
  REQUIRE(mat.searchBlockUID(1, 0) == 1);
  REQUIRE(mat.searchBlockUID(2, 0) == 2);

  REQUIRE(mat.searchBlockUID(0, 1) == -1);
  REQUIRE(mat.searchBlockUID(1, 1) == -1);
  REQUIRE(mat.searchBlockUID(2, 1) == -1);

  REQUIRE(mat.searchBlockUID(0, 2) == -1);
  REQUIRE(mat.searchBlockUID(1, 2) == 3);
  REQUIRE(mat.searchBlockUID(2, 2) == 4);

  REQUIRE(mat.searchBlockUID(0, 3) == 5);
  REQUIRE(mat.searchBlockUID(1, 3) == 6);
  REQUIRE(mat.searchBlockUID(2, 3) == -1);

  // inefficient!
  for (int r = 0; r < mat.numBlocksRow(); r++) {
    for (int c = 0; c < mat.numBlocksCol(); c++) {
      if (mat.hasBlock(r, c)) {
        mat.block(r, c).setConstant(r + c / 1.0);
      }
    }
  }

  for (int r = 0; r < mat.numBlocksRow(); r++) {
    for (int c = 0; c < mat.numBlocksCol(); c++) {
      if (mat.hasBlock(r, c)) {
        REQUIRE((mat.block(r, c).array() == r + c / 1.0).all());
      }
    }
  }

  {
    auto it = mat.colBegin(0);
    REQUIRE(it() == 0);
    REQUIRE(it.row() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 1);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 2);
    REQUIRE(it.row() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = mat.colBegin(1);
    REQUIRE(it() == it.end());

    it = mat.colBegin(2);
    REQUIRE(it() == 3);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 4);
    REQUIRE(it.row() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = mat.colBegin(3);
    REQUIRE(it() == 5);
    REQUIRE(it.row() == 0);
    REQUIRE(it() != it.end());

    //auto & bij = it.block();
    //bij.setConstant(it.row() + 3 / 1.0);
    auto m = it.block();
    REQUIRE((it.block().array() == it.row() + 3 / 1.0).all());
    REQUIRE((m.array() == it.row() + 3 / 1.0).all());
    auto vold = it.block()(0, 0);
    it.block().setConstant(23);
    REQUIRE((it.block().array() == 23).all());
    it.block().setConstant(vold);
    REQUIRE((it.block().array() == it.row() + 3 / 1.0).all());

    it++;
    REQUIRE(it() == 6);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());
  }

  {
    const auto& matC = mat;
    auto it = matC.colBegin(0);
    REQUIRE(it() == 0);
    REQUIRE(it.row() == 0);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 1);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 2);
    REQUIRE(it.row() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = matC.colBegin(1);
    REQUIRE(it() == it.end());

    it = matC.colBegin(2);
    REQUIRE(it() == 3);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == 4);
    REQUIRE(it.row() == 2);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());

    it = matC.colBegin(3);
    REQUIRE(it() == 5);
    REQUIRE(it.row() == 0);
    REQUIRE(it() != it.end());

    auto m = it.block();
    REQUIRE((it.block().array() == it.row() + 3 / 1.0).all());
    // it.block().setConstant(23); // does not compile!! ok!

    it++;
    REQUIRE(it() == 6);
    REQUIRE(it.row() == 1);
    REQUIRE(it() != it.end());
    it++;
    REQUIRE(it() == it.end());
  }

}