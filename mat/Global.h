#pragma once

#include <cassert>

#define __MAT_ASSERT_FALSE() assert(false && "mat::how do you end up here?!?");

namespace mat {

  constexpr int Fixed = 0;
  constexpr int Dynamic = -1;
  constexpr int Variable = -2;

  constexpr int BlockDense = 1;
  constexpr int BlockDiagonal = 2;
  constexpr int BlockSparse = 3;
  constexpr int BlockCoeffSparse = 4;

  constexpr int ColMajor = 0;
  constexpr int RowMajor = 1;

  template<int Ordering>
  struct IsColMajor {
    static constexpr bool value = false;
  };
  template<>
  struct IsColMajor<ColMajor> {
    static constexpr bool value = true;
  };

  template<int Ordering>
  struct IsRowMajor {
    static constexpr bool value = false;
  };
  template<>
  struct IsRowMajor<RowMajor> {
    static constexpr bool value = true;
  };

}