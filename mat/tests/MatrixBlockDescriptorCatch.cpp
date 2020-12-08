#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/MatrixBlockDescriptor.hpp"


// let use W matrix in BA case
static constexpr int camSize = 6;
static constexpr int pointSize = 3;
static constexpr int numCams = 4;
static constexpr int numPoints = 10;

static const std::vector<int> camSizeVar = { 5,8,6,5 };
static const std::vector<int> pointSizeVar = { 2,2,3,3,4,4,5,5,1,1 };

TEMPLATE_TEST_CASE_SIG("MatrixBlockDescriptor", "[MatrixBlockDescriptor]", ((int BR, int BC, int NBR, int NBC), BR, BC, NBR, NBC),
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

  using MatDescT = mat::MatrixBlockDescriptor<BR, BC, NBR, NBC>;

  std::unique_ptr<MatDescT> mat;
  SECTION("default ctor") {
    mat.reset(new MatDescT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));
    mat->resize(camSize, numCams, pointSize, numPoints);

  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSize, pointSize));
    mat->resize(camSize, numCams, pointSize, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSize, numCams, pointSize, numPoints));
    mat->resize(camSize, numCams, pointSize, numPoints);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams);
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSize);
  }
}


TEMPLATE_TEST_CASE_SIG("MatrixBlockDescriptorVar", "[MatrixBlockDescriptor]", ((int NBR, int NBC), NBR, NBC),
  (numCams, numPoints),
  (numCams, mat::Dynamic),
  (mat::Dynamic, numPoints),
  (mat::Dynamic, mat::Dynamic)
)
{
  using MatDescT = mat::MatrixBlockDescriptor<mat::Variable, mat::Variable, NBR, NBC>;
  std::unique_ptr<MatDescT> mat;
  SECTION("default ctor") {
    mat.reset(new MatDescT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));
    mat->resize(camSizeVar, numCams, pointSizeVar, numPoints);

  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSizeVar, pointSizeVar));
    mat->resize(camSizeVar, numCams, pointSizeVar, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSizeVar, numCams, pointSizeVar, numPoints));
    mat->resize(camSizeVar, numCams, pointSizeVar, numPoints);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams); // vector are same size on purpose!
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSizeVar[i]);
  }
}
