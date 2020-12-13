#pragma once

#include "DimensionDescriptor.h"

#include "DimensionDescriptorTraits.hpp"

#include <cassert>

namespace mat {

  template<int NB>
  DimensionDescriptorBase<NB>::DimensionDescriptorBase() {

  }

  template<int NB>
  DimensionDescriptorBase<NB>::DimensionDescriptorBase(int nb)  {
    assert(nb == NB);
  }

  template<int NB>
  DimensionDescriptorBase<NB>::~DimensionDescriptorBase() {

  }

  template<int NB>
  void DimensionDescriptorBase<NB>::setNumBlocks(int nb) {
    assert(nb == NB);
  }

  //------------------------------------------------------------------------------------------------


  template<int B, int NB>
  DimensionDescriptor<B, NB>::DimensionDescriptor() {
    static_assert(B > 0, "DimensionDescriptor has to be templatized with B > 0 or B=mat::Dynamic (-1) or B= mat::Variable (-2)");
  }

  template<int B, int NB>
  DimensionDescriptor<B, NB>::DimensionDescriptor(int b) {
    static_assert(B > 0, "DimensionDescriptor has to be templatized with B > 0 or B=mat::Dynamic (-1) or B= mat::Variable (-2)");
    assert(b == B && "block size parameters differs from compile time, use BlockDescriptor<mat::Dynamic> instead");
  }

  template<int B, int NB>
  DimensionDescriptor<B, NB>::DimensionDescriptor(int b, int nBlocks)
    : DimensionDescriptorBase<NB>(nBlocks) {
    static_assert(B > 0, "DimensionDescriptor has to be templatized with B > 0 or B=mat::Dynamic (-1) or B= mat::Variable (-2)");
    assert(b == B && "block size parameters differs from compile time, use BlockDescriptor<mat::Dynamic> instead");
  }

  template<int B, int NB>
  DimensionDescriptor<B, NB>::~DimensionDescriptor() {

  }

  //------------------------------------------------------------------------------------------------

  template<int NB>
  DimensionDescriptor<mat::Dynamic, NB>::DimensionDescriptor()
    : _b(0)
  {

  }

  template<int NB>
  DimensionDescriptor<mat::Dynamic, NB>::DimensionDescriptor(int b)
    : _b(b)
  {
    assert(b > 0);
  }

  template<int NB>
  DimensionDescriptor<mat::Dynamic, NB>::DimensionDescriptor(int b, int nBlocks)
    : DimensionDescriptorBase<NB> (nBlocks), _b(b)
  {
    assert(b > 0);
  }

  template<int NB>
  DimensionDescriptor<mat::Dynamic, NB>::~DimensionDescriptor() {

  }

  //------------------------------------------------------------------------------------------------

  template<int NB>
  DimensionDescriptor<mat::Variable, NB>::DimensionDescriptor()
  {
    _bi.resize(this->numBlocks(), 0); // create N blocks, but empty.... no sense
    updateStarts();
  }

  template<int NB>
  DimensionDescriptor<mat::Variable, NB>::DimensionDescriptor(const std::vector<int>& bi)
    : DimensionDescriptorBase<NB>(int(bi.size())), _bi(bi)
  {
    for (int i = 0; i < bi.size(); i++) {
      assert(bi[i] > 0);
    }
    updateStarts();
  }

  template<int NB>
  DimensionDescriptor<mat::Variable, NB>::DimensionDescriptor(const std::vector<int>& bi, int nBlocks)
    : DimensionDescriptorBase<NB>(int(bi.size())), _bi(bi) {
    assert(int(bi.size()) == nBlocks);
    for (int i = 0; i < bi.size(); i++) {
      assert(bi[i] > 0);
    }
    updateStarts();
  }

  template<int NB>
  DimensionDescriptor<mat::Variable, NB>::~DimensionDescriptor() {

  }

  template<int NB>
  void DimensionDescriptor<mat::Variable, NB>::updateStarts() {
    _bi_start.resize(_bi.size() + 1);
    _bi_start[0] = 0;
    for (int i = 1; i < _bi_start.size(); i++) {
      _bi_start[i] = _bi_start[i - 1] + _bi[i - 1];
    }
  }

  template<int NB>
  void DimensionDescriptor<mat::Variable, NB>::resize(const std::vector<int>& bi, int nb)
  {
    assert(int(bi.size()) == nb);
    this->resize(bi);
  }

  template<int NB>
  void DimensionDescriptor<mat::Variable, NB>::resize(const std::vector<int>& bi) {
    this->_bi = bi;
    this->setNumBlocks(int(bi.size()));
    this->updateStarts();
  }

  template<int NB>
  void DimensionDescriptor<mat::Variable, NB>::addBlock(int b) {
    this->_bi.push_back(b);
    this->_bi_start.push_back(this->_bi_start.back() + b);
    this->setNumBlocks(int(_bi.size()));
  }

}