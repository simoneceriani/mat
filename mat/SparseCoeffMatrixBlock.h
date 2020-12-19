#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"
#include "SparsityPattern.h"
#include "SparsityPatternBlockDescriptor.h"

#include <type_traits>

namespace mat {


  template<class T, int Ordering, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class SparseCoeffMatrixBlock final : public MatrixBlockBase<BR, BC, NBR, NBC> {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockCoeffSparse, T, BR, BC, NBR, NBC>;

    using StorageType = typename Traits::template StorageType<Ordering>;
    using SubMatrixType = typename Traits::SubMatrixType;
    using BlockType = typename Traits::template BlockType<Ordering>;
    using ConstBlockType = typename Traits::template ConstBlockType<Ordering>;

    using BlockDescriptor = MatrixBlockDescriptor<BR, BC, NBR, NBC>;
    using DimensionDescriptorRow = typename BlockDescriptor::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename BlockDescriptor::DimensionDescriptorCol;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

  private:
    // 
    StorageType _mat;
    std::vector<int> _outerStarts;  // size = outersize
    std::vector<int> _innerIndexes; // size = nnz blocks
    std::vector<int> _uid2outer;  // size = nnz blocks, return the outer given uid

    // keep track on how to map raw sparse data to blocks
    std::shared_ptr<const SparsityPatternBlockDescriptor<Ordering>> _sparseCoeffMap;

    template <class OuterDesc, class InnerDesc>
    void populateSparseMat(const SparsityPattern<Ordering>& sp, const OuterDesc& outDesc, const InnerDesc& inDesc);

    void resizeImpl(const SparsityPattern<Ordering>& sp);

    // -1 if does not exist, otherwise num of element in the inner (not UID)
    int searchBlockInner(int r, int c) const;

    inline int outer(int r, int c) const {
      if (Ordering == mat::ColMajor) {
        return c;
      }
      else if (Ordering == mat::RowMajor) {
        return r;
      }
      else {
        ASSERT_FALSE();
      }
    }

    inline int row(int outer, int inner) const {
      if (Ordering == mat::ColMajor) {
        return inner;
      }
      else if (Ordering == mat::RowMajor) {
        return outer;
      }
      else {
        ASSERT_FALSE();
      }
    }

    inline int col(int outer, int inner) const {
      if (Ordering == mat::ColMajor) {
        return outer;
      }
      else if (Ordering == mat::RowMajor) {
        return inner;
      }
      else {
        ASSERT_FALSE();
      }
    }

    inline BlockType blockOuterInner(int o, int in, int r, int c) {
      return BlockType(_mat.coeffs().data() + _sparseCoeffMap->offset(o, in), this->rowBlockSize(r), this->colBlockSize(c), Eigen::OuterStride<>(_sparseCoeffMap->stride(o)));
    }

    inline ConstBlockType blockOuterInner(int o, int in, int r, int c) const {
      return ConstBlockType(_mat.coeffs().data() + _sparseCoeffMap->offset(o, in), this->rowBlockSize(r), this->colBlockSize(c), Eigen::OuterStride<>(_sparseCoeffMap->stride(o)));
    }

  public:
    SparseCoeffMatrixBlock();
    SparseCoeffMatrixBlock(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering> &sp);

    virtual ~SparseCoeffMatrixBlock();

    void resize(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp);

    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    int nonZeroBlocks() const {
      return int(_innerIndexes.size());
    }

    int blockUID(int r, int c) const;

    bool hasBlock(int r, int c) {
      return searchBlockInner(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      int out = _uid2outer[uid];
      int in = uid - _outerStarts[out];
      int r = this->row(out, in);
      int c = this->col(out, in);
      return this->blockOuterInner(out, in, r, c);
    }

    inline ConstBlockType blockByUID(int uid) const {
      int out = _uid2outer[uid];
      int in = uid - _outerStarts[out];
      int r = this->row(out, in);
      int c = this->col(out, in);
      return this->blockOuterInner(out, in, r, c);
    }

    inline BlockType block(int r, int c) {
      int out = outer(r, c);
      int in = searchBlockInner(r, c);
      assert(in >= 0);
      return blockOuterInner(out, in, r, c);
    }

    inline ConstBlockType block(int r, int c) const {
      int out = outer(r, c);
      int in = searchBlockInner(r, c);
      assert(in >= 0);
      return blockOuterInner(out, in, r, c);
    }

    void setZero();

    // iterator on inner dimension 
    template<class BaseT>
    class InnerIterator {
      int _outer;

      BaseT* _sm;
      int _curId;
      int _lastId;

    public:
      InnerIterator(BaseT & _sm, int outer, int curId, int lastId);
      virtual ~InnerIterator();
      
      inline int operator()() const {
        assert(_curId <= _lastId);
        return _curId;
      }

      inline int blockUID() const {
        assert(_curId <= _lastId);
        return _curId;
      }

      inline int end() const {
        return _lastId;
      }

      inline void operator ++() {
        assert(_curId <= _lastId);
        _curId++;
      }
      inline void operator ++(int /*trash*/) {
        assert(_curId <= _lastId);
        _curId++;
      }

      
      int row() const {
        assert(_curId <= _lastId);
        int inner = _sm->_innerIndexes[_curId];
        return _sm->row(_outer, inner);
      }

      
      int col() const {
        assert(_curId <= _lastId);
        int inner = _sm->_innerIndexes[_curId];
        return _sm->col(_outer, inner);
      }

      template<typename RetType = BlockType>
      inline std::enable_if_t<!std::is_const<BaseT>::value, RetType> block() {
        assert(_curId <= _lastId);
        int inner = _sm->_innerIndexes[_curId];
        int r = _sm->row(_outer, inner);
        int c = _sm->col(_outer, inner);
        return _sm->blockOuterInner(_outer, inner, r, c);
      }

      inline ConstBlockType block() const {
        assert(_curId <= _lastId);
        int inner = _sm->_innerIndexes[_curId];
        int r = _sm->row(_outer, inner);
        int c = _sm->col(_outer, inner);
        return _sm->blockOuterInner(_outer, inner, r, c);
      }


    };

    // iterate on column blocks (if block are stored in col major)
    template<typename RetType = InnerIterator<SparseCoeffMatrixBlock>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) {
      return InnerIterator<SparseCoeffMatrixBlock>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    template<typename RetType = const InnerIterator<const SparseCoeffMatrixBlock>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) const {
      return InnerIterator<const SparseCoeffMatrixBlock>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    // iterate on row blocks (if block are stored in row major)
    template<typename RetType = InnerIterator<SparseCoeffMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) {
      return InnerIterator<SparseCoeffMatrixBlock>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }

    template<typename RetType = InnerIterator<const SparseCoeffMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) const {
      return InnerIterator<const SparseCoeffMatrixBlock>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }

  };

}