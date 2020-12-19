#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"
#include "SparsityPattern.h"

#include <type_traits>

namespace mat {

  template<class T, int Ordering, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class SparseMatrixBlock final : public MatrixBlockBase<BR, BC, NBR, NBC> {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockSparse, T, BR, BC, NBR, NBC>;

    using StorageType = typename Traits::StorageType;
    using SubMatrixType = typename Traits::SubMatrixType;
    using BlockType = typename Traits::BlockType;
    using ConstBlockType = typename Traits::ConstBlockType;

    using BlockDescriptor = MatrixBlockDescriptor<BR, BC, NBR, NBC>;
    using DimensionDescriptorRow = typename BlockDescriptor::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename BlockDescriptor::DimensionDescriptorCol;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

  private:
    // store all the blocks
    StorageType _mat;
    std::vector<int> _outerStarts;
    std::vector<int> _innerIndexes;

    void resizeImpl(const SparsityPattern<Ordering>& sp);

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
  public:
    SparseMatrixBlock();
    SparseMatrixBlock(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp);

    virtual ~SparseMatrixBlock();

    void resize(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp);

    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    int nonZeroBlocks() const {
      return int(_mat.size());
    }

    // -1 if does not exist
    int blockUID(int r, int c) const;

    bool hasBlock(int r, int c) const {
      return blockUID(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      return _mat[uid];
    }

    inline ConstBlockType blockByUID(int uid) const {
      return _mat[uid];
    }

    inline BlockType block(int r, int c) {
      return _mat[blockUID(r, c)];
    }

    inline ConstBlockType block(int r, int c) const {
      return _mat[blockUID(r, c)];
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
      InnerIterator(BaseT& _sm, int outer, int curId, int lastId);
      virtual ~InnerIterator();

      inline int operator()() const {
        assert(_curId <= _lastId);
        return _curId;
      }

      inline int blockUID() const {
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
        return _sm->_mat[_curId];
      }

      inline ConstBlockType block() const {
        assert(_curId <= _lastId);
        return _sm->_mat[_curId];
      }


    };

    // iterate on column blocks (if block are stored in col major)
    template<typename RetType = InnerIterator<SparseMatrixBlock>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) {
      return InnerIterator<SparseMatrixBlock>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    template<typename RetType = const InnerIterator<const SparseMatrixBlock>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) const {
      return InnerIterator<const SparseMatrixBlock>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    // iterate on row blocks (if block are stored in row major)
    template<typename RetType = InnerIterator<SparseMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) {
      return InnerIterator<SparseMatrixBlock>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }

    template<typename RetType = InnerIterator<const SparseMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) const {
      return InnerIterator<const SparseMatrixBlock>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }

  };

}