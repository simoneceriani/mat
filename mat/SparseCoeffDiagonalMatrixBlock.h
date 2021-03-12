#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "MatrixBlockBase.h"
#include "SparsityPattern.h"


namespace mat {

  template<class T, int Ordering, int BR, int BC, int NBR = mat::Dynamic, int NBC = mat::Dynamic >
  class SparseCoeffDiagonalMatrixBlock final : public MatrixBlockBase<BR, BC, NBR, NBC> {

  public:

    using Traits = MatrixBlockTypeTraits<mat::SparseCoeffBlockDiagonal, T, BR, BC, NBR, NBC>;

    using StorageType = typename Traits::template StorageType<Ordering>;
    using SubMatrixType = typename Traits::SubMatrixType;
    using BlockType = typename Traits::template BlockType<Ordering>;;
    using ConstBlockType = typename Traits::template ConstBlockType<Ordering>;

    using BlockDescriptor = MatrixBlockDescriptor<BR, BC, NBR, NBC>;
    using DimensionDescriptorRow = typename BlockDescriptor::DimensionDescriptorRow;
    using DimensionDescriptorCol = typename BlockDescriptor::DimensionDescriptorCol;

    using RowTraits = typename DimensionDescriptorRow::Traits;
    using ColTraits = typename DimensionDescriptorCol::Traits;

  private:
    StorageType _mat;
    Eigen::VectorXi _coeffOffsets;  // size = nnz blocks

    inline BlockType blockImpl(int i) {
      return BlockType(_mat.coeffs().data() + _coeffOffsets(i), this->rowBlockSize(i), this->colBlockSize(i));
    }

    inline ConstBlockType blockImpl(int i) const {
      return ConstBlockType(_mat.coeffs().data() + _coeffOffsets(i), this->rowBlockSize(i), this->colBlockSize(i));
    }

    typename SparsityPattern<Ordering>::CSPtr _sparsityPattern;

    template <class OuterDesc, class InnerDesc>
    void populateSparseMat(const OuterDesc& outDesc, const InnerDesc& inDesc, int nnz);

    template <class OuterDesc, class InnerDesc>
    int computeStrides(const OuterDesc& outDesc, const InnerDesc& inDesc);

    void resizeImpl();

  public:
    SparseCoeffDiagonalMatrixBlock();
    SparseCoeffDiagonalMatrixBlock(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp);

    virtual ~SparseCoeffDiagonalMatrixBlock();

    void resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp);

    const StorageType& mat() const {
      return _mat;
    }

    StorageType& mat() {
      return _mat;
    }

    int nonZeroBlocks() const {
      return int(_coeffOffsets.size());
    }

    int blockUID(int r, int c) const {
      if (r == c && r < this->nonZeroBlocks()) {
        return r;
      }
      return -1;
    }

    bool hasBlock(int r, int c) {
      return blockUID(r, c) >= 0;
    }

    inline BlockType blockByUID(int uid) {
      return this->blockImpl(uid);
    }

    inline ConstBlockType blockByUID(int uid) const {
      return this->blockImpl(uid);
    }

    inline BlockType block(int r, int c) {
      return this->blockImpl(r);
    }

    inline ConstBlockType block(int r, int c) const {
      return this->blockImpl(r);
    }

    void setZero();

    const typename SparsityPattern<Ordering>::CSPtr& sparsityPatternCSPtr() const {
      return _sparsityPattern;
    }

    const SparsityPattern<Ordering>& sparsityPattern() const {
      return *_sparsityPattern;
    }

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
    template<class BaseT>
    class InnerIterator {
      int _curId;
      BaseT* _sm;

      int _lastId;

    public:
      InnerIterator(BaseT& sm, int id, int lastId);
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

      int outerSize() const {
        return nonZeroBlocks();
      }

      // iterate on inner, be careful if it is rows or cols depends on ordering (no matter in this case, diagonal)
      InnerIterator<SparseCoeffDiagonalMatrixBlock> begin(int o) {
        return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
      }

      const InnerIterator<const SparseCoeffDiagonalMatrixBlock> begin(int o) const {
        return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
      }


      InnerIterator<SparseCoeffDiagonalMatrixBlock> colBegin(int c) {
        return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, c, (c < this->nonZeroBlocks() ? c + 1 : c));
      }

      InnerIterator<const SparseCoeffDiagonalMatrixBlock> colBegin(int c) const {
        return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, c, (c < this->nonZeroBlocks() ? c + 1 : c));
      }

      InnerIterator<SparseCoeffDiagonalMatrixBlock> rowBegin(int c) {
        return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, c, (c < this->nonZeroBlocks() ? c + 1 : c));
      }

      InnerIterator<const SparseCoeffDiagonalMatrixBlock> rowBegin(int c) const {
        return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, c, (c < this->nonZeroBlocks() ? c + 1 : c));
      }


    };

    // iterate on inner, be careful if it is rows or cols depends on ordering (no matter in this case, diagonal)    
    InnerIterator<SparseCoeffDiagonalMatrixBlock> begin(int o) {
      return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }
    const InnerIterator<const SparseCoeffDiagonalMatrixBlock> begin(int o) const {
      return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }

    // iterate on column blocks (if block are stored in col major)
    
    InnerIterator<SparseCoeffDiagonalMatrixBlock> colBegin(int o) {
      return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }

    InnerIterator<SparseCoeffDiagonalMatrixBlock> colBegin(int o) const {
      return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }

    // iterate on row blocks (if block are stored in row major)
    InnerIterator<SparseCoeffDiagonalMatrixBlock> rowBegin(int o) {
      return InnerIterator<SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }

    InnerIterator<SparseCoeffDiagonalMatrixBlock> rowBegin(int o) const {
      return InnerIterator<const SparseCoeffDiagonalMatrixBlock>(*this, o, (o < this->nonZeroBlocks() ? o + 1 : o));
    }

    inline int blockCoeffStart(int outer, int inner) const {
      return _coeffOffsets(outer);
    }
    int blockCoeffStride(int outer) const;

  };

}