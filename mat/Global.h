#pragma once

namespace mat {

  constexpr int Fixed = 0;
  constexpr int Dynamic = -1;
  constexpr int Variable = -2;

  constexpr int BlockDense = 1;
  constexpr int BlockDiagonal = 2;
  constexpr int BlockSparse = 3;

}