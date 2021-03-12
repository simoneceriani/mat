#pragma once
#include "SparseCoeffDiagonalMatrixBlock.h"

#include "MatrixBlockDescriptor.hpp"
#include "MatrixBlockBase.hpp"

namespace mat {

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseCoeffDiagonalMatrixBlock() {
    auto sp = std::make_shared<SparsityPattern<Ordering>>(this->numBlocksRow(), this->numBlocksCol());
    sp->setDiagonal();
    this->_sparsityPattern = sp;
    if (this->numBlocksRow() != 0 && this->numBlocksCol() != 0) {
      this->resizeImpl();
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::SparseCoeffDiagonalMatrixBlock(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp)
    : MatrixBlockBase<BR, BC, NBR, NBC>(blockDesc),
    _sparsityPattern(sp)
  {
#ifndef NDEBUG
    // check sparse pattern is diagonal (debug only)
    for (int o = 0; o < std::min(sp->outerSize(), sp->innerSize()); o++) {
      const auto& ins = sp->inner(o);
      assert(ins.size() == 1);
      assert(*(ins.begin()) == o);
    }
#endif
    if (blockDesc.numBlocksRow() != 0 && blockDesc.numBlocksCol() != 0) {
      this->resizeImpl();
    }
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::~SparseCoeffDiagonalMatrixBlock() {

  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template <class OuterDesc, class InnerDesc>
  int SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::computeStrides(const OuterDesc& outDesc, const InnerDesc& inDesc) {
    _coeffOffsets.resize(std::min(outDesc.numBlocks(), inDesc.numBlocks()));

    int offset = 0;
    if (InnerDesc::Traits::blockType == mat::Variable || OuterDesc::Traits::blockType == mat::Variable) {
      for (int oi = 0; oi < _coeffOffsets.size(); oi++) {
        _coeffOffsets(oi) = offset;
        offset += (inDesc.blockSize(oi) * outDesc.blockSize(oi));
      }
    }
    else {
      const int bs = inDesc.uniqueBlockSize() * outDesc.uniqueBlockSize();
      offset = _coeffOffsets.size() * bs;
      for (int oi = 0; oi < _coeffOffsets.size(); oi++) {
        _coeffOffsets(oi) = oi * bs;
      }
    }
    return offset;
  }


  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template <class OuterDesc, class InnerDesc>
  void SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::populateSparseMat(const OuterDesc& outDesc, const InnerDesc& inDesc, int nnz) {
    _mat.resize(this->numElementRows(), this->numElementCols());
    _mat.reserve(nnz);

    // Comments are relative to colmajor
    for (int o = 0; o < std::min(outDesc.numBlocks(), inDesc.numBlocks()); o++) {
      // o is the col index
      const int o_start = outDesc.blockStart(o);
      // iterate on the colums of this block, oi is the number of col in the block
      for (int oi = 0; oi < outDesc.blockSize(o); oi++) {
        _mat.startVec(o_start + oi);

        // the block is the diagonal one
        const int in = o;
        const int in_start = inDesc.blockStart(in);
        // iterate on rows of this block
        for (int ini = 0; ini < inDesc.blockSize(in); ini++) {
          // add element to sparse matrix
          _mat.insertBackByOuterInner(o_start + oi, in_start + ini) = 0; // init with 0
        }

      }
    }

    _mat.finalize();
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resizeImpl() {

    const auto& bd = this->blockDescriptor();
    // populate the matrix
    if (Ordering == mat::ColMajor) {
      int nnz = computeStrides(bd.colDescription(), bd.rowDescription());
      populateSparseMat(this->blockDescriptor().colDescription(), this->blockDescriptor().rowDescription(), nnz);
    }
    else if (Ordering == mat::RowMajor) {
      int nnz = computeStrides(bd.rowDescription(), bd.colDescription());
      populateSparseMat(this->blockDescriptor().rowDescription(), this->blockDescriptor().colDescription(), nnz);
    }
    else __MAT_ASSERT_FALSE();

  }


  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::resize(const BlockDescriptor& blockDesc, const typename SparsityPattern<Ordering>::CSPtr& sp) {
    MatrixBlockBase<BR, BC, NBR, NBC>::resize(blockDesc);
    _sparsityPattern = sp;
    this->resizeImpl();
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  void SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::setZero() {
    _mat.coeffs().setZero();
  }

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  int SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::blockCoeffStride(int outer) const {
    if (Ordering == mat::ColMajor) {
      if (ColTraits::blockType == mat::Variable || RowTraits::blockType == mat::Variable) {
        const auto& bd = this->blockDescriptor();
        return bd.colDescription().blockSize(outer);
      }
      else {
        const auto& bd = this->blockDescriptor();
        return bd.colDescription().uniqueBlockSize();
      }
    }
    else if (Ordering == mat::RowMajor) {
      if (ColTraits::blockType == mat::Variable || RowTraits::blockType == mat::Variable) {
        const auto& bd = this->blockDescriptor();
        return bd.rowDescription().blockSize(outer);
      }
      else {
        const auto& bd = this->blockDescriptor();
        return bd.rowDescription().uniqueBlockSize();
      }
    }
    else __MAT_ASSERT_FALSE();


  }

  //----------------------------------------------------------------------------

  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::InnerIterator(BaseT& sm, int id, int lastId) :
    _curId(id), _sm(&sm), _lastId(lastId) {

  }
  template< class T, int Ordering, int BR, int BC, int NBR, int NBC >
  template<class BaseT>
  SparseCoeffDiagonalMatrixBlock<T, Ordering, BR, BC, NBR, NBC>::InnerIterator<BaseT>::~InnerIterator() {

  }

}