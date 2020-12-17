#include "SparsityPatternBlockDescriptor.h"
#include "MatrixBlockDescriptor.hpp"
#include "SparsityPattern.h"

namespace mat {


  template<int Ordering>
  template<int BR, int BC, int NBR, int NBC>
  SparsityPatternBlockDescriptor<Ordering>::SparsityPatternBlockDescriptor(const SparsityPattern<Ordering>& sp, const MatrixBlockDescriptor<BR, BC, NBR, NBC>& bd)
    : _outerStrides(sp.outerSize()),
    _offsets(sp.outerSize()),
    _nonZeroCoeffs(0)
  {

    if (Ordering == mat::ColMajor) {
      computeStrides(sp, bd.colDescription(), bd.rowDescription());
    }
    else if (Ordering == mat::RowMajor) {
      computeStrides(sp, bd.rowDescription(), bd.colDescription());
    }
    else ASSERT_FALSE();

  }

  template<int Ordering>
  SparsityPatternBlockDescriptor< Ordering>::~SparsityPatternBlockDescriptor() {

  }

  template<int Ordering>
  template <class OuterDesc, class InnerDesc>
  void SparsityPatternBlockDescriptor< Ordering>::computeStrides(const SparsityPattern<Ordering>& sp, const OuterDesc& outDesc, const InnerDesc& inDesc) {
    int offset = 0;
    for (int oi = 0; oi < _outerStrides.size(); oi++) {
      const auto& blocksInner = sp.inner(oi);

      if (InnerDesc::Traits::blockType == mat::Variable) {
        int coeff_stride = 0;
        int offset_first = offset;
        int num_el = 0;
        int count = 0;
        _offsets[oi].resize(blocksInner.size());
        for (auto bi : blocksInner) {
          _offsets[oi][count++] = offset;
          int bs = inDesc.blockSize(bi);
          offset += bs;
          coeff_stride += bs;

          num_el += inDesc.blockSize(bi) * outDesc.blockSize(oi);
        }
        _outerStrides[oi] = coeff_stride;
        // prepare offset for next outer block
        offset = offset_first + num_el;
        // cumulate num elements
        _nonZeroCoeffs += num_el;
      }
      else {
        assert(inDesc.uniqueBlockSize() >= 0);
        _outerStrides[oi] = blocksInner.size() * inDesc.uniqueBlockSize();

        int offset_first = offset;
        _offsets[oi].resize(blocksInner.size());
        int count = 0;
        for (auto bi : blocksInner) {
          _offsets[oi][count++] = offset;
          offset += inDesc.uniqueBlockSize();
        }
        int num_el = blocksInner.size() * inDesc.uniqueBlockSize() * outDesc.blockSize(oi);
        offset = offset_first + num_el;
        // cumulate num elements
        _nonZeroCoeffs += num_el;
      }
    }

  }


}