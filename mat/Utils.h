#pragma once

namespace mat {

  class Utils {

  public:

    template <class ForwardIterator>
    static ForwardIterator binary_search(ForwardIterator first, ForwardIterator last, int val);

  };
}