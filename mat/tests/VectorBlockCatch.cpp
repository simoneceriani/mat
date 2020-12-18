#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/VectorBlock.hpp"

const int nr = 10;

static constexpr int pointSize = 3;
static constexpr int numPoints = 10;

static const std::vector<int> pointSizeVar = { 2,2,3,3,4,4,5,5,1,1 };

TEMPLATE_TEST_CASE_SIG("VectorBlock", "[VectorBlock]", ((int BR, int NBR), BR, NBR),
  (pointSize, numPoints),
  (pointSize, mat::Dynamic),
  (mat::Dynamic, numPoints),
  (mat::Dynamic, mat::Dynamic)
)
{

  using MatT = mat::VectorBlock<double, BR, NBR>;

  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numSegments() == (NBR == mat::Dynamic ? 0 : NBR));

    REQUIRE(mat->mat().rows() == (NBR == mat::Dynamic || BR == mat::Dynamic ? 0 : pointSize * numPoints));
    REQUIRE(mat->mat().cols() == 1);

    mat->resize(pointSize, numPoints);

  }  
  SECTION("sized ctor") {
    mat.reset(new MatT(pointSize));
    mat->resize(pointSize, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(pointSize, numPoints));
    mat->resize(pointSize, numPoints);
  }

  REQUIRE(mat->numElements() == pointSize * numPoints);
  REQUIRE(1 == mat->mat().cols());
  REQUIRE(mat->numElements() == mat->mat().rows());

  MatT mat2(mat->segmentDescriptionCSPtr());
  REQUIRE(mat2.numElements() == mat->mat().rows());
  REQUIRE(1== mat->mat().cols());
  REQUIRE(mat2.numElements() == mat2.mat().rows());
  REQUIRE(1 == mat2.mat().cols());

  for (int i = 0; i < mat->numSegments(); i++) {
    REQUIRE(mat->segmentStart(i) == i * pointSize);
    REQUIRE(mat->segmentSize(i) == pointSize);
  }
}

TEMPLATE_TEST_CASE_SIG("VectorBlockVar", "[VectorBlock]", ((int NBR), NBR),
  (numPoints),
  (mat::Dynamic)
)
{
  using MatT = mat::VectorBlock<double, mat::Variable, NBR>;
  std::unique_ptr<MatT> mat;
  SECTION("default ctor") {
    mat.reset(new MatT());
    REQUIRE(mat->numSegments() == (NBR == mat::Dynamic ? 0 : NBR));

    REQUIRE(mat->mat().rows() == 0);
    REQUIRE(mat->mat().cols() == 1);

    mat->resize(pointSizeVar, numPoints);
  }

  SECTION("sized ctor") {
    mat.reset(new MatT(pointSizeVar));
    mat->resize(pointSizeVar, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatT(pointSizeVar, numPoints));
    mat->resize(pointSizeVar, numPoints);
  }

  REQUIRE(mat->numElements() == pointSize * numPoints);
  REQUIRE(1 == mat->mat().cols());
  REQUIRE(mat->numElements() == mat->mat().rows());

  MatT mat2(mat->segmentDescriptionCSPtr());
  REQUIRE(mat2.numElements() == mat->mat().rows());
  REQUIRE(1 == mat->mat().cols());
  REQUIRE(mat2.numElements() == mat2.mat().rows());
  REQUIRE(1 == mat2.mat().cols());

  for (int i = 0; i < mat->numSegments(); i++) {
    REQUIRE(mat->segmentSize(i) == pointSizeVar[i]);
  }
}
