#pragma once
#include "Global.h"

#include <Eigen/Core>
#include <vector>

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

    void setNumBlocks(int nb);

  };

  //-----------------------------------------------------------------------------------------

  template<int B, int NB>
  struct DimensionDescriptorTraits {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = B;
    static constexpr int numElementsAtCompileTime = NB * B;
  };

  template<int B>
  struct DimensionDescriptorTraits<B, mat::Dynamic> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = B;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<int NB>
  struct DimensionDescriptorTraits<mat::Dynamic, NB> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = mat::Dynamic;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<>
  struct DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic> {
    using BlockSizeType = int;
    using BlockSizeTypePar = int;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = mat::Dynamic;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<int NB>
  struct DimensionDescriptorTraits<mat::Variable, NB> {
    using BlockSizeType = std::vector<int>;
    using BlockSizeTypePar = const std::vector<int>&;

    static constexpr int numBlocksAtCompileTime = NB;
    static constexpr int blockSizeAtCompileTime = mat::Variable;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
  };

  template<>
  struct DimensionDescriptorTraits<mat::Variable, mat::Dynamic> {
    using BlockSizeType = std::vector<int>;
    using BlockSizeTypePar = const std::vector<int>&;

    static constexpr int numBlocksAtCompileTime = mat::Dynamic;
    static constexpr int blockSizeAtCompileTime = mat::Variable;
    static constexpr int numElementsAtCompileTime = mat::Dynamic;
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

    inline int blockSize(int i) const {
      assert(i >= 0 && i < numBlocks());
      return B;
    }

    inline int blockStart(int i) const {
      assert(i >= 0 && i < numBlocks());
      return i * B;
    }

    void setBlockSize(int b) {
      assert(b == B);
    }

    void resize(int b, int nb) {
      this->setBlockSize(b);
      this->setNumBlocks(nb);
    }

    int numElements() const {
      return B * this->numBlocks();
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

    inline int blockSize(int i) const {
      assert(i >= 0 && i < numBlocks());
      return _b;
    }
    
    inline int blockStart(int i) const {
      assert(i >= 0 && i < numBlocks());
      return i * _b;
    }

    void setBlockSize(int b) {
      assert(b > 0);
      this->_b = b;
    }

    void resize(int b, int nb) {
      this->setBlockSize(b);
      this->setNumBlocks(nb);
    }

    int numElements() const {
      return this->_b * this->numBlocks();
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

    inline int blockSize(int i) const {
      assert(i >= 0 && i < numBlocks());
      return _bi[i];
    }

    inline int blockStart(int i) const {
      assert(i >= 0 && i < numBlocks());
      return _bi_start[i];
    }


    const std::vector<int>& blockSizes() const {
      return _bi;
    }

    void setBlockSize(const std::vector<int>& bi);
    void resize(const std::vector<int>& bi, int nb);

    void setNumBlocks(int nb);

    int numElements() const {
      return _bi_start.back();
    }

  };

}