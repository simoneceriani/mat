#pragma once
#include "Global.h"

#include <vector>

namespace mat {

  template<int Ordering> class SparsityPattern;
  template<int BR, int BC, int NBR, int NBC> class MatrixBlockDescriptor;

  template<int Ordering>
  class SparsityPatternBlockDescriptor {
    int _nonZeroCoeffs;
    std::vector<int> _outerStrides;
    std::vector<std::vector<int>> _offsets;

    template <class OuterDesc, class InnerDesc>
    void computeStrides(const SparsityPattern<Ordering>& sp, const OuterDesc& outDesc, const InnerDesc& inDesc);

  public:

    template<int BR, int BC, int NBR, int NBC>
    SparsityPatternBlockDescriptor(const SparsityPattern<Ordering>& sp, const MatrixBlockDescriptor<BR, BC, NBR, NBC>& bd);

    virtual ~SparsityPatternBlockDescriptor();

    int nonZeroCoeffs() const {
      return _nonZeroCoeffs;
    }

    int outerSize() const {
      return _outerStrides.size();
    }

    int stride(int outer) const {
      return _outerStrides[outer];
    }

    const std::vector<int>& outerStrides() const {
      return _outerStrides;
    }

    const std::vector<int>& offsets(int oi) const {
      return _offsets[oi];
    }

    int offset(int outer, int inner) const {
      return _offsets[outer][inner];
    }

  };

}