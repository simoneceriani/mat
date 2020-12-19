#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"
#include "SparsityPattern.h"


namespace mat {

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DiagonalMatrixBlock : public MatrixBlockBase<BR, BC, NBR, NBC>  {

  public:

    using Traits = MatrixBlockTypeTraits<mat::BlockDiagonal, T, BR, BC, NBR, NBC>;

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

    void resizeImpl();

  public:
    DiagonalMatrixBlock();

    // with block
    DiagonalMatrixBlock(const BlockDescriptor& blockDesc);

    virtual ~DiagonalMatrixBlock();

    void resize(const BlockDescriptor& blockDesc) override;


    const StorageType & mat() const {
      return _mat;
    }

    StorageType & mat() {
      return _mat;
    }

    inline BlockType block(int r, int c) {
      assert(r == c);
      assert(r >= 0 && r < _mat.size());
      return _mat[r];
    }

    inline ConstBlockType block(int r, int c) const {
      assert(r == c);
      assert(r >= 0 && r < _mat.size());
      return _mat[r];
    }

    void setZero();

  };

  template<class T, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class DiagonalMatrixBlockIterable final : public DiagonalMatrixBlock<T, BR, BC, NBR, NBC> {
  public:

    using DiagonalMatrixBlockT = DiagonalMatrixBlock<T, BR, BC, NBR, NBC>;

    using Traits = typename DiagonalMatrixBlockT::Traits;

    using StorageType = typename DiagonalMatrixBlockT::StorageType;
    using SubMatrixType = typename DiagonalMatrixBlockT::SubMatrixType;
    using BlockType = typename DiagonalMatrixBlockT::BlockType;
    using ConstBlockType = typename DiagonalMatrixBlockT::ConstBlockType;

    using BlockDescriptor = typename DiagonalMatrixBlockT::BlockDescriptor;
    using DimensionDescriptorRow = typename DiagonalMatrixBlockT::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename DiagonalMatrixBlockT::DimensionDescriptorCol;

    using RowTraits = typename DiagonalMatrixBlockT::RowTraits;
    using ColTraits = typename DiagonalMatrixBlockT::ColTraits;

  public:
    DiagonalMatrixBlockIterable();

    // with block
    template<int Ordering>
    DiagonalMatrixBlockIterable(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp);

    virtual ~DiagonalMatrixBlockIterable();

    template<int Ordering>
    void resize(const BlockDescriptor& blockDesc, const SparsityPattern<Ordering>& sp);

    int nonZeroBlocks() const {
      return this->mat().size();
    }

    int blockUID(int r, int c) const {
      if (r == c) {
        return r;
      }
      return -1;
    }

    bool hasBlock(int r, int c) {
      return blockUID(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      return this->block(uid,uid);
    }

    inline ConstBlockType blockByUID(int uid) const {
      return this->block(uid, uid);
    }

    template<class BaseT>
    class InnerIterator {
      int _curId;
      BaseT* _sm;

      int _lastId;

    public:
      InnerIterator(BaseT& sm, int id);
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
        return _curId;
      }

      int col() const {
        assert(_curId <= _lastId);
        return _curId;
      }

      template<typename RetType = BlockType>
      inline std::enable_if_t<!std::is_const<BaseT>::value, RetType> block() {
        assert(_curId <= _lastId);
        return _sm->blockByUID(_curId);
      }

      inline ConstBlockType block() const {
        assert(_curId <= _lastId);
        return _sm->blockByUID(_curId);
      }


    };

    InnerIterator<DiagonalMatrixBlockIterable> colBegin(int c) {
      return InnerIterator<DiagonalMatrixBlockIterable>(*this, c);
    }

    InnerIterator<const DiagonalMatrixBlockIterable> colBegin(int c) const {
      return InnerIterator<const DiagonalMatrixBlockIterable>(*this, c);
    }

    InnerIterator<DiagonalMatrixBlockIterable> rowBegin(int c) {
      return InnerIterator<DiagonalMatrixBlockIterable>(*this, c);
    }

    InnerIterator<const DiagonalMatrixBlockIterable> rowBegin(int c) const {
      return InnerIterator<const DiagonalMatrixBlockIterable>(*this, c);
    }

  };
}