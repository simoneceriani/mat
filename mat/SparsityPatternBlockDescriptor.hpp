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
    using MBD = MatrixBlockDescriptor<BR, BC, NBR, NBC>;


    int offset = 0;
    for (int oi = 0; oi < _outerStrides.size(); oi++) {
      const auto& blocksInner = sp.inner(oi);
      if (Ordering == mat::ColMajor) {

        if (MBD::RowTraits::blockType == mat::Variable) {
          int coeff_stride = 0;
          int offset_first = offset;
          int num_el = 0;
          int count = 0;
          _offsets[oi].resize(blocksInner.size());
          for (auto bi : blocksInner) {
            _offsets[oi][count++] = offset;
            int bs = bd.rowBlockSize(bi);
            offset += bs;
            coeff_stride += bs;

            num_el += bd.rowBlockSize(bi) * bd.colBlockSize(oi);
          }
          _outerStrides[oi] = coeff_stride;
          // prepare offset for next outer block
          offset = offset_first + num_el;
          // cumulate num elements
          _nonZeroCoeffs += num_el;
        }
        else {
          assert(bd.rowUniqueBlockSize() >= 0);
          _outerStrides[oi] = blocksInner.size() * bd.rowUniqueBlockSize();

          int offset_first = offset;
          _offsets[oi].resize(blocksInner.size());
          int count = 0;
          for (auto bi : blocksInner) {
            _offsets[oi][count++] = offset;
            offset += bd.rowUniqueBlockSize();
          }
          int num_el = blocksInner.size() * bd.rowUniqueBlockSize() * bd.colBlockSize(oi);
          offset = offset_first + num_el;
          // cumulate num elements
          _nonZeroCoeffs += num_el;
        }
      }
      else if (Ordering == mat::RowMajor) {
        if (MBD::ColTraits::blockType == mat::Variable) {
          int coeff_stride = 0;
          int offset_first = offset;
          int num_el = 0;
          int count = 0;
          _offsets[oi].resize(blocksInner.size());
          for (auto bi : blocksInner) {
            _offsets[oi][count++] = offset;
            int bs = bd.colBlockSize(bi);
            offset += bs;
            coeff_stride += bs;

            num_el += bd.colBlockSize(bi) * bd.rowBlockSize(oi);
          }
          _outerStrides[oi] = coeff_stride;
          // prepare offset for next outer block
          offset = offset_first + num_el;
          // cumulate num elements
          _nonZeroCoeffs += num_el;
        }
        else {
          assert(bd.colUniqueBlockSize() >= 0);
          _outerStrides[oi] = blocksInner.size() * bd.colUniqueBlockSize();

          int offset_first = offset;
          _offsets[oi].resize(blocksInner.size());
          int count = 0;
          for (auto bi : blocksInner) {
            _offsets[oi][count++] = offset;
            offset += bd.colUniqueBlockSize();
          }
          int num_el = blocksInner.size() * bd.colUniqueBlockSize() * bd.rowBlockSize(oi);
          offset = offset_first + num_el;
          // cumulate num elements
          _nonZeroCoeffs += num_el;
        }
      }
      else {
        ASSERT_FALSE();
      }
    }
  }

  template<int Ordering>
  SparsityPatternBlockDescriptor< Ordering>::~SparsityPatternBlockDescriptor() {

  }

}