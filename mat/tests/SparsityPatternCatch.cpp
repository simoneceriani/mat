#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/SparsityPattern.h"

const int nr = 10;
const int nc = 13;

TEST_CASE("SparsityPattern-RowMajor", "[SparsityPattern]") {
  mat::SparsityPattern<mat::RowMajor> sp(nr, nc);
  int count = 0;
  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
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

  for (int i = 0; i < nr; i++) {
    const auto& cols = sp.inner(i);
    for (auto c : cols) {
      REQUIRE(sp.has(i, c));
      REQUIRE(((i + c) % 2) == 0);
    }
  }

  auto mat = sp.toSparseMatrix();
  REQUIRE(mat.nonZeros() == sp.count());
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      if (sp.has(i, j)) {
        REQUIRE(mat.coeff(i, j) == 1);
      }
      else {
        REQUIRE(mat.coeff(i, j) == 0);
      }
    }
  }

  REQUIRE(mat.outerSize() == sp.outerSize());
  REQUIRE(mat.innerSize() == sp.innerSize());

  for (int k = 0; k < mat.outerSize(); ++k) {
    auto its = sp.inner(k).begin();
    for (Eigen::SparseMatrix<int, Eigen::RowMajor>::InnerIterator it(mat, k); it; ++it)
    {
      REQUIRE(it.value() == 1);
      REQUIRE(it.index() == *its);
      its++;
    }
  }
}

TEST_CASE("SparsityPattern-ColMajor", "[SparsityPattern]") {
  mat::SparsityPattern<mat::ColMajor> sp(nr, nc);
  int count = 0;
  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nc; j++) {
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
  for (int j = 0; j < nc; j++) {
    const auto& rows = sp.inner(j);
    for (auto r : rows) {
      REQUIRE(sp.has(r, j));
      REQUIRE(((r + j) % 2) == 0);
    }
  }

  auto mat = sp.toSparseMatrix();
  REQUIRE(mat.nonZeros() == sp.count());
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      if (sp.has(i, j)) {
        REQUIRE(mat.coeff(i, j) == 1);
      }
      else {
        REQUIRE(mat.coeff(i, j) == 0);
      }
    }
  }

  REQUIRE(mat.outerSize() == sp.outerSize());
  REQUIRE(mat.innerSize() == sp.innerSize());
  for (int k = 0; k < mat.outerSize(); ++k) {
    auto its = sp.inner(k).begin();
    for (Eigen::SparseMatrix<int, Eigen::ColMajor>::InnerIterator it(mat, k); it; ++it)
    {
      REQUIRE(it.value() == 1);
      REQUIRE(it.index() == *its);
      its++;
    }
  }

}
