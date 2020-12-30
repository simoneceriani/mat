#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"
#include "SparsityPattern.h"

namespace mat {

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DenseMatrixBlock : public MatrixBlockBase<BR, BC, NBR, NBC> {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockDense, T, BR, BC, NBR, NBC>;

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
    StorageType _mat;


  public:
    DenseMatrixBlock();

    // with block
    DenseMatrixBlock(const BlockDescriptor& blockDesc);

    virtual ~DenseMatrixBlock();

    virtual void resize(const BlockDescriptor& blockDesc) override final;


    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    inline BlockType block(int r, int c) {
      return _mat.template block<RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>(
        this->rowBlockStart(r), this->colBlockStart(c),
        this->rowBlockSize(r), this->colBlockSize(c)
        );
    }

    inline ConstBlockType block(int r, int c) const {
      return _mat.template block<RowTraits::blockSizeAtCompileTime, ColTraits::blockSizeAtCompileTime>(
        this->rowBlockStart(r), this->colBlockStart(c),
        this->rowBlockSize(r), this->colBlockSize(c)
        );
    }

    void setZero();

  };

  template<class T, int Ordering, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DenseMatrixBlockIterable final : public DenseMatrixBlock<T, BR, BC, NBR, NBC> {

  public:

    using DenseMatrixBlockT = DenseMatrixBlock<T, BR, BC, NBR, NBC>;

    using Traits = typename DenseMatrixBlockT::Traits;

    using StorageType = typename DenseMatrixBlockT::StorageType;
    using SubMatrixType = typename DenseMatrixBlockT::SubMatrixType;
    using BlockType = typename DenseMatrixBlockT::BlockType;
    using ConstBlockType = typename DenseMatrixBlockT::ConstBlockType;

    using BlockDescriptor = typename DenseMatrixBlockT::BlockDescriptor;
    using DimensionDescriptorRow = typename DenseMatrixBlockT::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename DenseMatrixBlockT::DimensionDescriptorCol;

    using RowTraits = typename DenseMatrixBlockT::RowTraits;
    using ColTraits = typename DenseMatrixBlockT::ColTraits;

  private:
    std::vector<int> _outerStarts;  // size = outersize
    std::vector<int> _innerIndexes; // size = nnz blocks
    std::vector<int> _uid2outer;  // size = nnz blocks, return the outer given uid

    typename SparsityPattern<Ordering>::CSPtr _sparsityPattern;

    void createPattern(const SparsityPattern<Ordering>& sp);
    
    inline int row(int outer, int inner) const {
      if (Ordering == mat::ColMajor) {
        return inner;
      }
      else if (Ordering == mat::RowMajor) {
        return outer;
      }
      else {
        __MAT_ASSERT_FALSE();
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
        __MAT_ASSERT_FALSE();
      }
    }
  public:
    DenseMatrixBlockIterable();
    DenseMatrixBlockIterable(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp);

    virtual ~DenseMatrixBlockIterable();

    void resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp);

    int nonZeroBlocks() const {
      return int(_innerIndexes.size());
    }

    int blockUID(int r, int c) const;

    bool hasBlock(int r, int c) {
      return blockUID(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      int out = _uid2outer[uid];
      int in = uid - _outerStarts[out];
      int r = this->row(out, in);
      int c = this->col(out, in);
      return this->block(r, c);
    }

    inline ConstBlockType blockByUID(int uid) const {
      int out = _uid2outer[uid];
      int in = uid - _outerStarts[out];
      int r = this->row(out, in);
      int c = this->col(out, in);
      return this->block(r, c);
    }

    const typename SparsityPattern<Ordering>::CSPtr & sparsityPatternCSPtr() const {
      return _sparsityPattern;
    }

    const SparsityPattern<Ordering> & sparsityPattern() const {
      return * _sparsityPattern;
    }

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
        return _sm->block(r, c);
      }

      inline ConstBlockType block() const {
        assert(_curId <= _lastId);
        int inner = _sm->_innerIndexes[_curId];
        int r = _sm->row(_outer, inner);
        int c = _sm->col(_outer, inner);
        return _sm->block(r, c);
      }


    };

    // iterate on column blocks (if block are stored in col major)
    template<typename RetType = InnerIterator<DenseMatrixBlockIterable>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) {
      return InnerIterator<DenseMatrixBlockIterable>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    template<typename RetType = const InnerIterator<const DenseMatrixBlockIterable>>
    std::enable_if_t<IsColMajor<Ordering>::value, RetType> colBegin(int c) const {
      return InnerIterator<const DenseMatrixBlockIterable>(*this, c, _outerStarts[c], _outerStarts[c + 1]);
    }

    // iterate on row blocks (if block are stored in row major)
    template<typename RetType = InnerIterator<DenseMatrixBlockIterable>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) {
      return InnerIterator<DenseMatrixBlockIterable>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }

    template<typename RetType = InnerIterator<const DenseMatrixBlockIterable>>
    std::enable_if_t<IsRowMajor<Ordering>::value, RetType> rowBegin(int r) const {
      return InnerIterator<const DenseMatrixBlockIterable>(*this, r, _outerStarts[r], _outerStarts[r + 1]);
    }
  };
}