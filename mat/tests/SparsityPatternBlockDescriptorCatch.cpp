#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/SparsityPatternBlockDescriptor.hpp"


// let use W matrix in BA case
static constexpr int camSize = 6;
static constexpr int pointSize = 3;
static constexpr int numCams = 3;
static constexpr int numPoints = 5;

static const std::vector<int> camSizeVar = { 5,7,6 };
static const std::vector<int> pointSizeVar = { 2,3,4,5,1 };

TEMPLATE_TEST_CASE_SIG("SparsityPatternBlockDescriptor", "[SparsityPatternBlockDescriptor]", ((int BR, int BC, int NBR, int NBC, int Ordering), BR, BC, NBR, NBC, Ordering),
  (camSize, pointSize, numCams, numPoints, mat::ColMajor),
  (camSize, pointSize, numCams, mat::Dynamic, mat::ColMajor),
  (camSize, pointSize, mat::Dynamic, numPoints, mat::ColMajor),
  (camSize, pointSize, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (camSize, mat::Dynamic, numCams, numPoints, mat::ColMajor),
  (camSize, mat::Dynamic, numCams, mat::Dynamic, mat::ColMajor),
  (camSize, mat::Dynamic, mat::Dynamic, numPoints, mat::ColMajor),
  (camSize, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, pointSize, numCams, numPoints, mat::ColMajor),
  (mat::Dynamic, pointSize, numCams, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, pointSize, mat::Dynamic, numPoints, mat::ColMajor),
  (mat::Dynamic, pointSize, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, numCams, numPoints, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, numCams, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, numPoints, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (camSize, pointSize, numCams, numPoints, mat::RowMajor),
  (camSize, pointSize, numCams, mat::Dynamic, mat::RowMajor),
  (camSize, pointSize, mat::Dynamic, numPoints, mat::RowMajor),
  (camSize, pointSize, mat::Dynamic, mat::Dynamic, mat::RowMajor),
  (camSize, mat::Dynamic, numCams, numPoints, mat::RowMajor),
  (camSize, mat::Dynamic, numCams, mat::Dynamic, mat::RowMajor),
  (camSize, mat::Dynamic, mat::Dynamic, numPoints, mat::RowMajor),
  (camSize, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, pointSize, numCams, numPoints, mat::RowMajor),
  (mat::Dynamic, pointSize, numCams, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, pointSize, mat::Dynamic, numPoints, mat::RowMajor),
  (mat::Dynamic, pointSize, mat::Dynamic, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, numCams, numPoints, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, numCams, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, numPoints, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::RowMajor)
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

  mat::SparsityPattern<Ordering> sp(numCams, numPoints);
  /*
  +-------------------+
  | X |   | X | X |   |
  +-------------------+
  |   | X |   |   | X |
  +-------------------+
  | X | X |   |   | X |
  */
  sp.add(0, 0); sp.add(0, 2); sp.add(0, 3);
  sp.add(1, 1); sp.add(1, 4);
  sp.add(2, 0); sp.add(2, 1); sp.add(2, 4);

  mat::SparsityPatternBlockDescriptor<Ordering> spbd(sp, *mat);

  if (Ordering == mat::ColMajor) {
    REQUIRE(sp.outerSize() == numPoints);
    REQUIRE(spbd.outerSize() == numPoints);

    REQUIRE(spbd.stride(0) == mat->rowBlockSize(0) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->rowBlockSize(0));

    int numEl = mat->colBlockSize(0) * (mat->rowBlockSize(0) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(1) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(1) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(2) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(2, 0) == numEl);
    numEl += mat->colBlockSize(2) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(3) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(3, 0) == numEl);
    numEl += mat->colBlockSize(3) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(4) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(4, 0) == numEl);
    REQUIRE(spbd.offset(4, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(4) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(numEl == spbd.nonZeroCoeffs());

  }
  else if (Ordering == mat::RowMajor) {
    REQUIRE(sp.outerSize() == numCams);
    REQUIRE(spbd.outerSize() == numCams);

    REQUIRE(spbd.stride(0) == mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->colBlockSize(0));
    REQUIRE(spbd.offset(0, 2) == mat->colBlockSize(0) + mat->colBlockSize(2));

    int numEl = mat->rowBlockSize(0) * (mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.stride(1) == mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(1) * (mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(spbd.stride(2) == mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(2, 0) == numEl);
    REQUIRE(spbd.offset(2, 1) == numEl + mat->colBlockSize(0));
    REQUIRE(spbd.offset(2, 2) == numEl + mat->colBlockSize(0) + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(2) * (mat->colBlockSize(0)+ mat->colBlockSize(1)+ mat->colBlockSize(4));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }


}


TEMPLATE_TEST_CASE_SIG("SparsityPatternBlockDescriptor-Var", "[MatrixBlockDescriptor]", ((int NBR, int NBC, int Ordering), NBR, NBC, Ordering),
  (numCams, numPoints, mat::ColMajor),
  (numCams, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, numPoints, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (numCams, numPoints, mat::RowMajor),
  (numCams, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, numPoints, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, mat::RowMajor)
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

  mat::SparsityPattern<Ordering> sp(numCams, numPoints);
  /*
  +-------------------+
  | X |   | X | X |   |
  +-------------------+
  |   | X |   |   | X |
  +-------------------+
  | X | X |   |   | X |
  */
  sp.add(0, 0); sp.add(0, 2); sp.add(0, 3);
  sp.add(1, 1); sp.add(1, 4);
  sp.add(2, 0); sp.add(2, 1); sp.add(2, 4);

  mat::SparsityPatternBlockDescriptor<Ordering> spbd(sp, *mat);

  if (Ordering == mat::ColMajor) {
    REQUIRE(sp.outerSize() == numPoints);
    REQUIRE(spbd.outerSize() == numPoints);

    REQUIRE(spbd.stride(0) == mat->rowBlockSize(0) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->rowBlockSize(0));

    int numEl = mat->colBlockSize(0) * (mat->rowBlockSize(0) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(1) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(1) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(2) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(2, 0) == numEl);
    numEl += mat->colBlockSize(2) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(3) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(3, 0) == numEl);
    numEl += mat->colBlockSize(3) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(4) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(4, 0) == numEl);
    REQUIRE(spbd.offset(4, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(4) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
  else if (Ordering == mat::RowMajor) {
    REQUIRE(sp.outerSize() == numCams);
    REQUIRE(spbd.outerSize() == numCams);

    REQUIRE(spbd.stride(0) == mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->colBlockSize(0));
    REQUIRE(spbd.offset(0, 2) == mat->colBlockSize(0) + mat->colBlockSize(2));

    int numEl = mat->rowBlockSize(0) * (mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.stride(1) == mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(1) * (mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(spbd.stride(2) == mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(2, 0) == numEl);
    REQUIRE(spbd.offset(2, 1) == numEl + mat->colBlockSize(0));
    REQUIRE(spbd.offset(2, 2) == numEl + mat->colBlockSize(0) + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(2) * (mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
}

TEMPLATE_TEST_CASE_SIG("SparsityPatternBlockDescriptor-VarFix", "[MatrixBlockDescriptor]", ((int BC, int NBR, int NBC, int Ordering), BC, NBR, NBC, Ordering),
  (pointSize, numCams, numPoints, mat::ColMajor),
  (pointSize, numCams, mat::Dynamic, mat::ColMajor),
  (pointSize, mat::Dynamic, numPoints, mat::ColMajor),
  (pointSize, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (pointSize, numCams, numPoints, mat::RowMajor),
  (pointSize, numCams, mat::Dynamic, mat::RowMajor),
  (pointSize, mat::Dynamic, numPoints, mat::RowMajor),
  (pointSize, mat::Dynamic, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, numCams, numPoints, mat::ColMajor),
  (mat::Dynamic, numCams, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, numPoints, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, numCams, numPoints, mat::RowMajor),
  (mat::Dynamic, numCams, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, numPoints, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::RowMajor)
)
{
  using MatDescT = mat::MatrixBlockDescriptor<mat::Variable, BC, NBR, NBC>;
  std::unique_ptr<MatDescT> mat;
  SECTION("default ctor") {
    mat.reset(new MatDescT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));
    mat->resize(camSizeVar, numCams, pointSize, numPoints);

  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSizeVar, pointSize));
    mat->resize(camSizeVar, numCams, pointSize, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSizeVar, numCams, pointSize, numPoints));
    mat->resize(camSizeVar, numCams, pointSize, numPoints);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams); // vector are same size on purpose!
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSizeVar[i]);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSize);
  }

  mat::SparsityPattern<Ordering> sp(numCams, numPoints);
  /*
  +-------------------+
  | X |   | X | X |   |
  +-------------------+
  |   | X |   |   | X |
  +-------------------+
  | X | X |   |   | X |
  */
  sp.add(0, 0); sp.add(0, 2); sp.add(0, 3);
  sp.add(1, 1); sp.add(1, 4);
  sp.add(2, 0); sp.add(2, 1); sp.add(2, 4);

  mat::SparsityPatternBlockDescriptor<Ordering> spbd(sp, *mat);

  if (Ordering == mat::ColMajor) {
    REQUIRE(sp.outerSize() == numPoints);
    REQUIRE(spbd.outerSize() == numPoints);

    REQUIRE(spbd.stride(0) == mat->rowBlockSize(0) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->rowBlockSize(0));

    int numEl = mat->colBlockSize(0) * (mat->rowBlockSize(0) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(1) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(1) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(2) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(2, 0) == numEl);
    numEl += mat->colBlockSize(2) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(3) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(3, 0) == numEl);
    numEl += mat->colBlockSize(3) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(4) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(4, 0) == numEl);
    REQUIRE(spbd.offset(4, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(4) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
  else if (Ordering == mat::RowMajor) {
    REQUIRE(sp.outerSize() == numCams);
    REQUIRE(spbd.outerSize() == numCams);

    REQUIRE(spbd.stride(0) == mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->colBlockSize(0));
    REQUIRE(spbd.offset(0, 2) == mat->colBlockSize(0) + mat->colBlockSize(2));

    int numEl = mat->rowBlockSize(0) * (mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.stride(1) == mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(1) * (mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(spbd.stride(2) == mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(2, 0) == numEl);
    REQUIRE(spbd.offset(2, 1) == numEl + mat->colBlockSize(0));
    REQUIRE(spbd.offset(2, 2) == numEl + mat->colBlockSize(0) + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(2) * (mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
}

TEMPLATE_TEST_CASE_SIG("SparsityPatternBlockDescriptor-FixVar", "[MatrixBlockDescriptor]", ((int BR, int NBR, int NBC, int Ordering), BR, NBR, NBC, Ordering),
  (camSize, numCams, numPoints, mat::ColMajor),
  (camSize, numCams, mat::Dynamic, mat::ColMajor),
  (camSize, mat::Dynamic, numPoints, mat::ColMajor),
  (camSize, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (camSize, numCams, numPoints, mat::RowMajor),
  (camSize, numCams, mat::Dynamic, mat::RowMajor),
  (camSize, mat::Dynamic, numPoints, mat::RowMajor),
  (camSize, mat::Dynamic, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, numCams, numPoints, mat::ColMajor),
  (mat::Dynamic, numCams, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, numPoints, mat::ColMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::ColMajor),
  (mat::Dynamic, numCams, numPoints, mat::RowMajor),
  (mat::Dynamic, numCams, mat::Dynamic, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, numPoints, mat::RowMajor),
  (mat::Dynamic, mat::Dynamic, mat::Dynamic, mat::RowMajor)
)
{
  using MatDescT = mat::MatrixBlockDescriptor<BR, mat::Variable, NBR, NBC>;
  std::unique_ptr<MatDescT> mat;
  SECTION("default ctor") {
    mat.reset(new MatDescT());
    REQUIRE(mat->numBlocksRow() == (NBR == mat::Dynamic ? 0 : NBR));
    REQUIRE(mat->numBlocksCol() == (NBC == mat::Dynamic ? 0 : NBC));
    mat->resize(camSize, numCams, pointSizeVar, numPoints);

  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSize, pointSizeVar));
    mat->resize(camSize, numCams, pointSizeVar, numPoints);
  }
  SECTION("sized ctor") {
    mat.reset(new MatDescT(camSize, numCams, pointSizeVar, numPoints));
    mat->resize(camSize, numCams, pointSizeVar, numPoints);
  }

  REQUIRE(mat->numElementRows() == camSize * numCams); // vector are same size on purpose!
  REQUIRE(mat->numElementCols() == pointSize * numPoints);

  for (int i = 0; i < mat->numBlocksRow(); i++) {
    REQUIRE(mat->rowBlockSize(i) == camSize);
  }
  for (int i = 0; i < mat->numBlocksCol(); i++) {
    REQUIRE(mat->colBlockSize(i) == pointSizeVar[i]);
  }

  mat::SparsityPattern<Ordering> sp(numCams, numPoints);
  /*
  +-------------------+
  | X |   | X | X |   |
  +-------------------+
  |   | X |   |   | X |
  +-------------------+
  | X | X |   |   | X |
  */
  sp.add(0, 0); sp.add(0, 2); sp.add(0, 3);
  sp.add(1, 1); sp.add(1, 4);
  sp.add(2, 0); sp.add(2, 1); sp.add(2, 4);

  mat::SparsityPatternBlockDescriptor<Ordering> spbd(sp, *mat);

  if (Ordering == mat::ColMajor) {
    REQUIRE(sp.outerSize() == numPoints);
    REQUIRE(spbd.outerSize() == numPoints);

    REQUIRE(spbd.stride(0) == mat->rowBlockSize(0) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->rowBlockSize(0));

    int numEl = mat->colBlockSize(0) * (mat->rowBlockSize(0) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(1) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(1) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(spbd.stride(2) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(2, 0) == numEl);
    numEl += mat->colBlockSize(2) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(3) == mat->rowBlockSize(0));
    REQUIRE(spbd.offset(3, 0) == numEl);
    numEl += mat->colBlockSize(3) * mat->rowBlockSize(0);

    REQUIRE(spbd.stride(4) == mat->rowBlockSize(1) + mat->rowBlockSize(2));
    REQUIRE(spbd.offset(4, 0) == numEl);
    REQUIRE(spbd.offset(4, 1) == numEl + mat->rowBlockSize(1));
    numEl += mat->colBlockSize(4) * (mat->rowBlockSize(1) + mat->rowBlockSize(2));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
  else if (Ordering == mat::RowMajor) {
    REQUIRE(sp.outerSize() == numCams);
    REQUIRE(spbd.outerSize() == numCams);

    REQUIRE(spbd.stride(0) == mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.offset(0, 0) == 0);
    REQUIRE(spbd.offset(0, 1) == mat->colBlockSize(0));
    REQUIRE(spbd.offset(0, 2) == mat->colBlockSize(0) + mat->colBlockSize(2));

    int numEl = mat->rowBlockSize(0) * (mat->colBlockSize(0) + mat->colBlockSize(2) + mat->colBlockSize(3));

    REQUIRE(spbd.stride(1) == mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(1, 0) == numEl);
    REQUIRE(spbd.offset(1, 1) == numEl + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(1) * (mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(spbd.stride(2) == mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));
    REQUIRE(spbd.offset(2, 0) == numEl);
    REQUIRE(spbd.offset(2, 1) == numEl + mat->colBlockSize(0));
    REQUIRE(spbd.offset(2, 2) == numEl + mat->colBlockSize(0) + mat->colBlockSize(1));
    numEl += mat->rowBlockSize(2) * (mat->colBlockSize(0) + mat->colBlockSize(1) + mat->colBlockSize(4));

    REQUIRE(numEl == spbd.nonZeroCoeffs());
  }
}
