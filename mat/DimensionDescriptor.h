#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

#include "DimensionDescriptorTraits.h"

namespace mat {

  template<int NB>
  class DimensionDescriptorBase {

  public:
    DimensionDescriptorBase();
    DimensionDescriptorBase(int nb);
    virtual ~DimensionDescriptorBase();

    inline int numBlocks() const {
      return NB;
    }

  protected:
    void setNumBlocks(int nb);

  };
  
  //-----------------------------------------------------------------------------------------

  template<>
  class DimensionDescriptorBase<mat::Dynamic> {
    int _nb;

  public:
    DimensionDescriptorBase();
    DimensionDescriptorBase(int nb);
    virtual ~DimensionDescriptorBase();

    inline int numBlocks() const {
      return _nb;
    }

  protected:
    void setNumBlocks(int nb);

  };

  //-----------------------------------------------------------------------------------------

  template<int B, int NB>
  class DimensionDescriptor : public DimensionDescriptorBase<NB> {

  public:

    using Traits = DimensionDescriptorTraits<B, NB>;

    DimensionDescriptor();
    DimensionDescriptor(int b);
    DimensionDescriptor(int b, int nBlocks);
    virtual ~DimensionDescriptor();

    inline int uniqueBlockSize() const {
      return B;
    }

    inline int blockSize(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return B;
    }

    inline int blockStart(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return i * B;
    }

    void resize(int b, int nb) {
      assert(b == B);
      this->setNumBlocks(nb);
    }

    int numElements() const {
      return B * this->numBlocks();
    }

    void addBlock(int b) {
      assert(b == B);
      this->setNumBlocks(this->numBlocks() + 1);
    }

  };

  //-----------------------------------------------------------------------------------------

  template<int NB>
  class DimensionDescriptor<mat::Dynamic, NB> : public DimensionDescriptorBase<NB> {

    int _b;

  public:

    using Traits = DimensionDescriptorTraits<mat::Dynamic, NB>;

    DimensionDescriptor();
    DimensionDescriptor(int b);
    DimensionDescriptor(int b, int nBlocks);
    virtual ~DimensionDescriptor();

    inline int uniqueBlockSize() const {
      return _b;
    }

    inline int blockSize(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return _b;
    }
    
    inline int blockStart(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return i * _b;
    }

    void resize(int b, int nb) {
      assert(b > 0);
      this->_b = b;
      this->setNumBlocks(nb);
    }

    int numElements() const {
      return this->_b * this->numBlocks();
    }

    void addBlock(int b) {
      this->setNumBlocks(this->numBlocks() + 1);
    }

  };

  //-----------------------------------------------------------------------------------------

  template<int NB>
  class DimensionDescriptor<mat::Variable, NB> : public DimensionDescriptorBase<NB> {

    std::vector<int> _bi;
    std::vector<int> _bi_start;

    void updateStarts();

  public:

    using Traits = DimensionDescriptorTraits<mat::Variable, NB>;

    DimensionDescriptor();
    DimensionDescriptor(const std::vector<int>& bi);
    DimensionDescriptor(const std::vector<int>& bi, int nBlocks);
    virtual ~DimensionDescriptor();

    inline int uniqueBlockSize() const {
      return -1;
    }

    inline int blockSize(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return _bi[i];
    }

    inline int blockStart(int i) const {
      assert(i >= 0 && i < this->numBlocks());
      return _bi_start[i];
    }


    const std::vector<int>& blockSizes() const {
      return _bi;
    }

    void resize(const std::vector<int>& bi, int nb);
    void resize(const std::vector<int>& bi);

    void addBlock(int b);


    int numElements() const {
      return _bi_start.back();
    }

  };

}