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

  public:
    SparseMatrixBlock();
    SparseMatrixBlock(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering> &sp);

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

    template <class ForwardIterator>
    static ForwardIterator binary_search(ForwardIterator first, ForwardIterator last, int val);

    // -1 if does not exist
    int searchBlockUID(int r, int c);

    bool hasBlock(int r, int c) {
      return searchBlockUID(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      return _mat[uid];
    }

    inline ConstBlockType blockByUID(int uid) const {
      return _mat[uid];
    }

    inline BlockType block(int r, int c) {
      return _mat[searchBlockUID(r, c)];
    }

    inline ConstBlockType block(int r, int c) const {
      return _mat[searchBlockUID(r, c)];
    }

    void setZero();

    // iterator on inner dimension 
    template<class BaseT>
    class InnerIterator {
      BaseT* _sm;
      int _curId;
      int _lastId;

    public:
      InnerIterator(BaseT & _sm, int curId, int lastId);
      virtual ~InnerIterator();
      
      inline int operator()() {
        assert(_curId <= _lastId);
        return _curId;
      }

      inline int end() {
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

      template<typename RetType = int>
      std::enable_if_t<IsColMajor<Ordering>::value, RetType> row() {
        assert(_curId <= _lastId);
        return _sm->_innerIndexes[_curId];
      }

      template<typename RetType = int>
      std::enable_if_t<IsRowMajor<Ordering>::value, RetType> col() {
        assert(_curId <= _lastId);
        return _sm->_innerIndexes[_curId];
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
      return InnerIterator<SparseMatrixBlock>(*this, _outerStarts[c], _outerStarts[c + 1]);
    }

    template<typename RetType = const InnerIterator<const SparseMatrixBlock>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) const {
      return InnerIterator<const SparseMatrixBlock>(*this, _outerStarts[c], _outerStarts[c + 1]);
    }

    // iterate on row blocks (if block are stored in row major)
    template<typename RetType = InnerIterator<SparseMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int c) {
      return InnerIterator<SparseMatrixBlock>(*this, _outerStarts[c], _outerStarts[c + 1]);
    }

    template<typename RetType = InnerIterator<const SparseMatrixBlock>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int c) const {
      return InnerIterator<const SparseMatrixBlock>(*this, _outerStarts[c], _outerStarts[c + 1]);
    }

  };

}