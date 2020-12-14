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
    mutable int _count;

  public:
    SparsityPattern(int nr, int nc);
    virtual ~SparsityPattern();

    void clear();

    void add(int i, int j) {
      if (Ordering == mat::ColMajor) {
        _sp[j].insert(i);
        _count = -1;
      }
      else if (Ordering == mat::RowMajor) {
        _sp[i].insert(j);
        _count = -1;
      }
      else {
        assert(false && "how do you end up here?!?");
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
        assert(false && "how do you end up here?!?");
      }
    }

    int count() const {
      if (_count == -1) {
        _count = 0;
        for (int k = 0; k < _sp.size(); k++) {
          _count += _sp[k].size();
        }
      }
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



}