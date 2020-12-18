#pragma once
#include "Utils.h"

#include <algorithm>

namespace mat {

  template <class ForwardIterator>
  static ForwardIterator Utils::binary_search(ForwardIterator first, ForwardIterator last, int val) {
    first = std::lower_bound(first, last, val);
    if (first != last) {
      if (val == *first) return first;
      return last;
    }
    return last;
  }

}