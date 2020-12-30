#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/DiagonalMatrixBlock.hpp"

// let use W matrix in BA case
static constexpr int camSize = 6;
static constexpr int pointSize = 3;
static constexpr int numCams = 4;
static constexpr int numPoints = 10;

static const std::vector<int> camSizeVar = { 5,8,6,5 };
static const std::vector<int> pointSizeVar = { 2,2,3,3,4,4,5,5,1,1 };

TEMPLATE_TEST_CASE_SIG("DiagonalMatrixBlockIterable", "[DiagonalMatrixBlockIterable]", ((int BR, int BC, int NBR, int NBC), BR, BC, NBR, NBC),
  (camSize, pointSize, numCams, numPoints),
  (camSize, pointSize, numCams, mat::Dynamic),
  (camSize, pointSize, mat::Dynamic, numPoints),
  (camSize, pointSize, mat::Dynamic, mat::Dynamic),
  (camSize, mat::Dynamic, numCams, numPoints),
  (camSize, mat::Dynamic, numCams, mat::Dynamic),
  (camSize, mat::Dynamic, mat::Dynamic, numPoints),
  (camSize, mat::Dynamic, mat::Dynamic, mat::Dynamic),
  (mat::Dynamic, pointSize, numCams, numPoints),
  (mat::Dynamic, pointSize, numCams, mat::Dynamic),
  (mat::Dynamic, pointSize, mat::Dynamic, numPoints),
  (mat::Dynamic, pointSize, mat::Dynamic, mat::Dynamic),
  (mat::Dynamic, mat::Dynamic, numCams, numPoints),
  (mat::Dynamic, mat::Dynamic, numCams, mat::Dynamic),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, numPoints),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::Dynamic)
)
{

  using MatT = mat::DiagonalMatrixBlockIterable<double, mat::ColMajor, BR, BC, NBR, NBC>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));

    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));

  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, pointSize), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  MatT mat2(mat->blockDescriptor(), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat2.numElementRows() == camSize * numCams);
  REQUIRE(mat2.numElementCols() == pointSize * numPoints);

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat3.numElementRows() == camSize * numCams);
  REQUIRE(mat3.numElementCols() == pointSize * numPoints);

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat4.numElementRows() == camSize * numCams);
  REQUIRE(mat4.numElementCols() == pointSize * numPoints);

  MatT mat5;
  mat5.resize(mat->blockDescriptor(), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat5.numElementRows() == camSize * numCams);
  REQUIRE(mat5.numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
    REQUIRE(mat2.rowBlockSize(i) == camSize);
    REQUIRE(mat3.rowBlockSize(i) == camSize);
    REQUIRE(mat4.rowBlockSize(i) == camSize);
    REQUIRE(mat5.rowBlockSize(i) == camSize);
    auto it = mat->rowBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());

  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSize);
    REQUIRE(mat2.colBlockSize(i) == pointSize);
    REQUIRE(mat3.colBlockSize(i) == pointSize);
    REQUIRE(mat4.colBlockSize(i) == pointSize);
    REQUIRE(mat5.colBlockSize(i) == pointSize);
    auto it = mat->colBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());
  }

}

TEMPLATE_TEST_CASE_SIG("DiagonalMatrixBlockIterable-Square", "[DiagonalMatrixBlockIterable]", ((int BR, int NBR), BR, NBR),
  (camSize, numCams),
  (camSize, mat::Dynamic),
  (mat::Dynamic, numCams),
  (mat::Dynamic, mat::Dynamic)
)
{

  using MatT = mat::DiagonalMatrixBlockIterable<double, mat::ColMajor, BR, BR, NBR, NBR>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBR == mat::Dynamic ? 0 : NBR));

    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));

  }
  SECTION("sized ctor") {
    mat.reset(new MatT(MatT::BlockDescriptor::squareMatrix(camSize), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams)));
    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(MatT::BlockDescriptor::squareMatrix(camSize, numCams), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams)));
    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == camSize * numCams);

  MatT mat2(MatT::BlockDescriptor::squareMatrix(mat->blockDescriptor().rowDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(mat2.numElementRows() == camSize * numCams);
  REQUIRE(mat2.numElementCols() == camSize * numCams);

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(mat3.numElementRows() == camSize * numCams);
  REQUIRE(mat3.numElementCols() == camSize * numCams);

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(mat4.numElementRows() == camSize * numCams);
  REQUIRE(mat4.numElementCols() == camSize * numCams);

  MatT mat5;
  mat5.resize(mat->blockDescriptor(), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(mat5.numElementRows() == camSize * numCams);
  REQUIRE(mat5.numElementCols() == camSize * numCams);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
    REQUIRE(mat2.rowBlockSize(i) == camSize);
    REQUIRE(mat3.rowBlockSize(i) == camSize);
    REQUIRE(mat4.rowBlockSize(i) == camSize);
    REQUIRE(mat5.rowBlockSize(i) == camSize);
    auto it = mat->rowBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == camSize);
    REQUIRE(mat2.colBlockSize(i) == camSize);
    REQUIRE(mat3.colBlockSize(i) == camSize);
    REQUIRE(mat4.colBlockSize(i) == camSize);
    REQUIRE(mat5.colBlockSize(i) == camSize);
    auto it = mat->colBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());
  }
}



TEMPLATE_TEST_CASE_SIG("DiagonalMatrixBlockIterable-Var", "[DiagonalMatrixBlockIterable]", ((int NBR, int NBC), NBR, NBC),
  (numCams, numPoints),
  (numCams, mat::Dynamic),
  (mat::Dynamic, numPoints),
  (mat::Dynamic, mat::Dynamic)
)
{
  using MatT = mat::DiagonalMatrixBlockIterable<double, mat::ColMajor, mat::Variable, mat::Variable, NBR, NBC>;
  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));

    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));

  }

  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSizeVar, pointSizeVar), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  MatT mat2(typename MatT::BlockDescriptor(mat->blockDescriptor()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat2.numElementRows() == camSize * numCams);
  REQUIRE(mat2.numElementCols() == pointSize * numPoints);

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat3.numElementRows() == camSize * numCams);
  REQUIRE(mat3.numElementCols() == pointSize * numPoints);

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat4.numElementRows() == camSize * numCams);
  REQUIRE(mat4.numElementCols() == pointSize * numPoints);

  MatT mat5;
  mat5.resize(typename MatT::BlockDescriptor(mat->blockDescriptor()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numPoints));
  REQUIRE(mat5.numElementRows() == camSize * numCams);
  REQUIRE(mat5.numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.rowBlockSize(i) == camSizeVar[i]);
    auto it = mat->rowBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat2.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat3.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat4.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat5.colBlockSize(i) == pointSizeVar[i]);
    auto it = mat->colBegin(i);
    for (; it() != it.end(); it++) {
      REQUIRE(it() == it.row());
      REQUIRE(it() == it.col());
      REQUIRE(it() == i);
    }
    REQUIRE(it() == it.end());
  }
}

TEMPLATE_TEST_CASE_SIG("DiagonalMatrixBlockIterable-Square", "[DiagonalMatrixBlockIterable]", ((int NBR), NBR),
  (numCams),
  (mat::Dynamic)
)
{
  using MatT = mat::DiagonalMatrixBlockIterable<double, mat::ColMajor, mat::Variable, mat::Variable, NBR, NBR>;
  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBR == mat::Dynamic ? 0 : NBR));

    REQUIRE(mat->mat().size() == (NBR == mat::Dynamic ? 0 : NBR));

    auto block = typename MatT::BlockDescriptor(camSizeVar, numCams, camSizeVar, numCams);
    mat->resize(block, mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));

  }

  SECTION("sized ctor") {
    auto block = MatT::BlockDescriptor::squareMatrix(camSizeVar);
    mat.reset(new MatT(block, mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams)));
    block.resizeSquare(camSizeVar, numCams);
    mat->resize(block, mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  }
  SECTION("sized ctor") {
    auto block = MatT::BlockDescriptor::squareMatrix(camSizeVar);
    mat.reset(new MatT(block, mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams)));
    block.resizeSquare(camSizeVar, numCams);
    mat->resize(block, mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == camSize * numCams);

  REQUIRE(mat->mat().size() == numCams);

  MatT mat2(mat->blockDescriptor(), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(std::min(mat2.numBlocksRow(), mat2.numBlocksCol()) == mat->mat().size());

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(std::min(mat3.numBlocksRow(), mat3.numBlocksCol()) == mat->mat().size());

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(std::min(mat4.numBlocksRow(), mat4.numBlocksCol()) == mat->mat().size());

  MatT mat5;
  mat5.resize(typename MatT::BlockDescriptor(mat->blockDescriptor()), mat::SparsityPatternColMajor::makeDiagonal(numCams, numCams));
  REQUIRE(std::min(mat5.numBlocksRow(), mat5.numBlocksCol()) == mat->mat().size());

  for (int i = 0; i < std::min(mat->numBlocksRow(), mat->numBlocksCol()); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.rowBlockSize(i) == camSizeVar[i]);

    REQUIRE(mat->colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.colBlockSize(i) == camSizeVar[i]);
  }
}
