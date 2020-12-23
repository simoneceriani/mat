#pragma once
#include "Global.h"

#include <set>
#include <vector>
#include <cassert>

#include <memory>

namespace mat {

  template<int Ordering>
  class SparsityPattern {
    std::vector<std::set<int>> _sp;
    int _count;

  public:
    SparsityPattern(int nr, int nc);
    virtual ~SparsityPattern();

    static SparsityPattern makeDiagonal(int nr, int nc);

    void clear();

    void setDiagonal() {
      for (int i = 0; i < _sp.size(); i++) {
        auto ret = _sp[i].insert(i);
        if (ret.second) _count++;
      }
    }

    void add(int i, int j) {
      if (Ordering == mat::ColMajor) {
        auto ret = _sp[j].insert(i);
        if (ret.second) _count++;
      }
      else if (Ordering == mat::RowMajor) {
        auto ret = _sp[i].insert(j);
        if (ret.second) _count++;
      }
      else {
        __MAT_ASSERT_FALSE();
      }
    }

    bool has(int i, int j) const {
      if (Ordering == mat::ColMajor) {
        return _sp[j].find(i) != _sp[j].end();
      }
      else if (Ordering == mat::RowMajor) {
        return _sp[i].find(j) != _sp[i].end();
      }
      else {
        __MAT_ASSERT_FALSE();
      }
    }

    int count() const {
      return _count;
    }

    // if rowmajor, outer is "i", return is the list of occupied cols
    // if colmajor, outer is "j", return is the list of occupied rows
    const std::set<int>& inner(int outer) const {
      return _sp[outer];
    }

    int outerSize() const {
      return _sp.size();
    }

  };

  using SparsityPatternColMajor = SparsityPattern<mat::ColMajor>;
  using SparsityPatternRowMajor = SparsityPattern<mat::RowMajor>;

}