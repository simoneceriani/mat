#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/DenseMatrixBlock.hpp"

// let use W matrix in BA case
static constexpr int camSize = 6;
static constexpr int pointSize = 3;
static constexpr int numCams = 4;
static constexpr int numPoints = 10;

static const std::vector<int> camSizeVar = { 5,8,6,5 };
static const std::vector<int> pointSizeVar = { 2,2,3,3,4,4,5,5,1,1 };

TEMPLATE_TEST_CASE_SIG("DenseMatrixBlock", "[DenseMatrixBlock]", ((int BR, int BC, int NBR, int NBC), BR, BC, NBR, NBC),
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

  using MatT = mat::DenseMatrixBlock<double,BR, BC, NBR, NBC>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));

    REQUIRE(mat->mat().rows() == (NBR == mat::Dynamic || BR == mat::Dynamic ? 0 : camSize * numCams));
    REQUIRE(mat->mat().cols() == (NBC == mat::Dynamic || BC == mat::Dynamic ? 0 : pointSize * numPoints));

    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints));

  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, pointSize)));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSize, numCams, pointSize, numPoints));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  REQUIRE(mat->numElementRows() == mat->mat().rows());
  REQUIRE(mat->numElementCols() == mat->mat().cols());

  MatT mat2(mat->blockDescriptor());
  REQUIRE(mat2.numElementRows() == mat->mat().rows());
  REQUIRE(mat2.numElementCols() == mat->mat().cols());
  REQUIRE(mat2.numElementRows() == mat2.mat().rows());
  REQUIRE(mat2.numElementCols() == mat2.mat().cols());

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat3.numElementRows() == mat->mat().rows());
  REQUIRE(mat3.numElementCols() == mat->mat().cols());
  REQUIRE(mat3.numElementRows() == mat3.mat().rows());
  REQUIRE(mat3.numElementCols() == mat3.mat().cols());

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat4.numElementRows() == mat->mat().rows());
  REQUIRE(mat4.numElementCols() == mat->mat().cols());
  REQUIRE(mat4.numElementRows() == mat4.mat().rows());
  REQUIRE(mat4.numElementCols() == mat4.mat().cols());

  MatT mat5;
  mat5.resize(mat->blockDescriptor());
  REQUIRE(mat5.numElementRows() == mat->mat().rows());
  REQUIRE(mat5.numElementCols() == mat->mat().cols());
  REQUIRE(mat5.numElementRows() == mat5.mat().rows());
  REQUIRE(mat5.numElementCols() == mat5.mat().cols());

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
    REQUIRE(mat2.rowBlockSize(i) == camSize);
    REQUIRE(mat3.rowBlockSize(i) == camSize);
    REQUIRE(mat4.rowBlockSize(i) == camSize);
    REQUIRE(mat5.rowBlockSize(i) == camSize);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSize);
    REQUIRE(mat2.colBlockSize(i) == pointSize);
    REQUIRE(mat3.colBlockSize(i) == pointSize);
    REQUIRE(mat4.colBlockSize(i) == pointSize);
    REQUIRE(mat5.colBlockSize(i) == pointSize);
  }
}

TEMPLATE_TEST_CASE_SIG("DenseMatrixBlock-Square", "[DenseMatrixBlock]", ((int BR, int NBR), BR, NBR),
  (camSize, numCams),
  (camSize, mat::Dynamic),
  (mat::Dynamic, numCams),
  (mat::Dynamic, mat::Dynamic)
)
{

  using MatT = mat::DenseMatrixBlock<double, BR, BR, NBR, NBR>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBR == mat::Dynamic ? 0 : NBR));

    REQUIRE(mat->mat().rows() == (NBR == mat::Dynamic || BR == mat::Dynamic ? 0 : camSize * numCams));
    REQUIRE(mat->mat().cols() == (NBR == mat::Dynamic || BR == mat::Dynamic ? 0 : camSize * numCams));

    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams));

  }
  SECTION("sized ctor") {
    mat.reset(new MatT(MatT::BlockDescriptor::squareMatrix(camSize)));
    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(MatT::BlockDescriptor::squareMatrix(camSize, numCams)));
    mat->resize(MatT::BlockDescriptor::squareMatrix(camSize, numCams));
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == camSize * numCams);

  REQUIRE(mat->numElementRows() == mat->mat().rows());
  REQUIRE(mat->numElementCols() == mat->mat().cols());

  MatT mat2(MatT::BlockDescriptor::squareMatrix(mat->blockDescriptor().rowDescriptionCSPtr()));
  REQUIRE(mat2.numElementRows() == mat->mat().rows());
  REQUIRE(mat2.numElementCols() == mat->mat().cols());
  REQUIRE(mat2.numElementRows() == mat2.mat().rows());
  REQUIRE(mat2.numElementCols() == mat2.mat().cols());

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat3.numElementRows() == mat->mat().rows());
  REQUIRE(mat3.numElementCols() == mat->mat().cols());
  REQUIRE(mat3.numElementRows() == mat3.mat().rows());
  REQUIRE(mat3.numElementCols() == mat3.mat().cols());

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat4.numElementRows() == mat->mat().rows());
  REQUIRE(mat4.numElementCols() == mat->mat().cols());
  REQUIRE(mat4.numElementRows() == mat4.mat().rows());
  REQUIRE(mat4.numElementCols() == mat4.mat().cols());

  MatT mat5;
  mat5.resize(mat->blockDescriptor());
  REQUIRE(mat5.numElementRows() == mat->mat().rows());
  REQUIRE(mat5.numElementCols() == mat->mat().cols());
  REQUIRE(mat5.numElementRows() == mat5.mat().rows());
  REQUIRE(mat5.numElementCols() == mat5.mat().cols());

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
    REQUIRE(mat2.rowBlockSize(i) == camSize);
    REQUIRE(mat3.rowBlockSize(i) == camSize);
    REQUIRE(mat4.rowBlockSize(i) == camSize);
    REQUIRE(mat5.rowBlockSize(i) == camSize);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == camSize);
    REQUIRE(mat2.colBlockSize(i) == camSize);
    REQUIRE(mat3.colBlockSize(i) == camSize);
    REQUIRE(mat4.colBlockSize(i) == camSize);
    REQUIRE(mat5.colBlockSize(i) == camSize);
  }
}



TEMPLATE_TEST_CASE_SIG("DenseMatrixBlockVar", "[DenseMatrixBlock]", ((int NBR, int NBC), NBR, NBC),
  (numCams, numPoints),
  (numCams, mat::Dynamic),
  (mat::Dynamic, numPoints),
  (mat::Dynamic, mat::Dynamic)
)
{
  using MatT = mat::DenseMatrixBlock<double, mat::Variable, mat::Variable, NBR, NBC>;
  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));

    REQUIRE(mat->mat().rows() == 0);
    REQUIRE(mat->mat().cols() == 0);

    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints));

  }
  
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSizeVar, pointSizeVar)));
    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints));
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints)));
    mat->resize(typename MatT::BlockDescriptor(camSizeVar, numCams, pointSizeVar, numPoints));
  }
  
  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  REQUIRE(mat->numElementRows() == mat->mat().rows());
  REQUIRE(mat->numElementCols() == mat->mat().cols());

  MatT mat2(typename MatT::BlockDescriptor(mat->blockDescriptor()));
  REQUIRE(mat2.numElementRows() == mat->mat().rows());
  REQUIRE(mat2.numElementCols() == mat->mat().cols());
  REQUIRE(mat2.numElementRows() == mat2.mat().rows());
  REQUIRE(mat2.numElementCols() == mat2.mat().cols());

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat3.numElementRows() == mat->mat().rows());
  REQUIRE(mat3.numElementCols() == mat->mat().cols());
  REQUIRE(mat3.numElementRows() == mat3.mat().rows());
  REQUIRE(mat3.numElementCols() == mat3.mat().cols());

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat4.numElementRows() == mat->mat().rows());
  REQUIRE(mat4.numElementCols() == mat->mat().cols());
  REQUIRE(mat4.numElementRows() == mat4.mat().rows());
  REQUIRE(mat4.numElementCols() == mat4.mat().cols());

  MatT mat5;
  mat5.resize(typename MatT::BlockDescriptor(mat->blockDescriptor()));
  REQUIRE(mat5.numElementRows() == mat->mat().rows());
  REQUIRE(mat5.numElementCols() == mat->mat().cols());
  REQUIRE(mat5.numElementRows() == mat5.mat().rows());
  REQUIRE(mat5.numElementCols() == mat5.mat().cols());

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.rowBlockSize(i) == camSizeVar[i]);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat2.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat3.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat4.colBlockSize(i) == pointSizeVar[i]);
    REQUIRE(mat5.colBlockSize(i) == pointSizeVar[i]);
  }
}

TEMPLATE_TEST_CASE_SIG("DenseMatrixBlockVar-Square", "[DenseMatrixBlock]", ((int NBR), NBR),
  (numCams),
  (mat::Dynamic)
)
{
  using MatT = mat::DenseMatrixBlock<double, mat::Variable, mat::Variable, NBR, NBR>;
  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBR == mat::Dynamic ? 0 : NBR));

    REQUIRE(mat->mat().rows() == 0);
    REQUIRE(mat->mat().cols() == 0);

    auto block = typename MatT::BlockDescriptor(camSizeVar, numCams, camSizeVar, numCams);
    mat->resize(block);

  }

  SECTION("sized ctor") {
    auto block = MatT::BlockDescriptor::squareMatrix(camSizeVar);
    mat.reset(new MatT(block));
    block.resizeSquare(camSizeVar, numCams);
    mat->resize(block);
  }
  SECTION("sized ctor") {
    auto block = MatT::BlockDescriptor::squareMatrix(camSizeVar);
    mat.reset(new MatT(block));
    block.resizeSquare(camSizeVar, numCams);
    mat->resize(block);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == camSize * numCams);

  REQUIRE(mat->numElementRows() == mat->mat().rows());
  REQUIRE(mat->numElementCols() == mat->mat().cols());

  MatT mat2(mat->blockDescriptor());
  REQUIRE(mat2.numElementRows() == mat->mat().rows());
  REQUIRE(mat2.numElementCols() == mat->mat().cols());
  REQUIRE(mat2.numElementRows() == mat2.mat().rows());
  REQUIRE(mat2.numElementCols() == mat2.mat().cols());

  MatT mat3(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat3.numElementRows() == mat->mat().rows());
  REQUIRE(mat3.numElementCols() == mat->mat().cols());
  REQUIRE(mat3.numElementRows() == mat3.mat().rows());
  REQUIRE(mat3.numElementCols() == mat3.mat().cols());

  MatT mat4;
  mat4.resize(typename MatT::BlockDescriptor(mat->blockDescriptor().rowDescriptionCSPtr(), mat->blockDescriptor().colDescriptionCSPtr()));
  REQUIRE(mat4.numElementRows() == mat->mat().rows());
  REQUIRE(mat4.numElementCols() == mat->mat().cols());
  REQUIRE(mat4.numElementRows() == mat4.mat().rows());
  REQUIRE(mat4.numElementCols() == mat4.mat().cols());

  MatT mat5;
  mat5.resize(typename MatT::BlockDescriptor(mat->blockDescriptor()));
  REQUIRE(mat5.numElementRows() == mat->mat().rows());
  REQUIRE(mat5.numElementCols() == mat->mat().cols());
  REQUIRE(mat5.numElementRows() == mat5.mat().rows());
  REQUIRE(mat5.numElementCols() == mat5.mat().cols());

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.rowBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.rowBlockSize(i) == camSizeVar[i]);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat2.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat3.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat4.colBlockSize(i) == camSizeVar[i]);
    REQUIRE(mat5.colBlockSize(i) == camSizeVar[i]);
  }
}
