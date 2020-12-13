#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "mat/DimensionDescriptor.hpp"

#include <memory>
#include <type_traits>

TEMPLATE_TEST_CASE_SIG("DimensionDescriptor<3|Dynamic>", "[DimensionDescriptor]", ((int B, int NB), B, NB), (3, 5), (3, mat::Dynamic), (mat::Dynamic, 5), (mat::Dynamic, mat::Dynamic)) {


  std::unique_ptr<mat::DimensionDescriptor<B, NB>> d;
  SECTION("Default ctor") {

    d.reset(new mat::DimensionDescriptor<B, NB>());
    REQUIRE(mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == NB);

    d->resize(3, 5);
    REQUIRE(d->numBlocks() == 5);


  }
  SECTION("Block-size ctor") {

    d.reset(new mat::DimensionDescriptor<B, NB>(3));
    REQUIRE(mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == NB);

    d->resize(3, 5);
    REQUIRE(d->numBlocks() == 5);


  }
  SECTION("Block-size + block num ctor") {

    d.reset(new mat::DimensionDescriptor<B, NB>(3, 5));
    REQUIRE(mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == NB);


    REQUIRE(d->numBlocks() == 5);
  }

  int s = 0;
  for (int i = 0; i < d->numBlocks(); i++) {
    REQUIRE(d->blockSize(i) == 3);
    REQUIRE(d->blockStart(i) == 3 * i);
    s += d->blockSize(i);
  }

  REQUIRE(d->numElements() == s);

  if (mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == mat::Dynamic) {
    int b = (mat::DimensionDescriptor<B, NB>::Traits::blockSizeAtCompileTime == mat::Dynamic ? 3 : mat::DimensionDescriptor<B, NB>::Traits::blockSizeAtCompileTime);
    d->addBlock(b);


    s = 0;
    for (int i = 0; i < d->numBlocks(); i++) {
      REQUIRE(d->blockSize(i) == 3);
      REQUIRE(d->blockStart(i) == s);
      s += d->blockSize(i);
    }

    REQUIRE(d->numElements() == s);
  }


}

TEMPLATE_TEST_CASE_SIG("DimensionDescriptor<mat::Variable|5>", "[DimensionDescriptor]", (int NB, NB), 5, mat::Dynamic) {
  constexpr int B = mat::Variable;

  std::vector<int> blockSizes = { 3,2,4,4,2 };

  std::unique_ptr<mat::DimensionDescriptor<B, NB>> d;
  SECTION("Default ctor") {

    d.reset(new mat::DimensionDescriptor<B, NB>());
    REQUIRE(mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == NB);
    

    d->resize(blockSizes, blockSizes.size());
    REQUIRE(d->numBlocks() == blockSizes.size());

    d->resize(blockSizes);
    REQUIRE(d->numBlocks() == blockSizes.size());


  }
  SECTION("Block-size ctor") {

    SECTION("vector only") {
      d.reset(new mat::DimensionDescriptor<B, NB>(blockSizes));
    }
    SECTION("vector + size only") {
      d.reset(new mat::DimensionDescriptor<B, NB>(blockSizes, int(blockSizes.size())));
    }

    REQUIRE(d->numBlocks() == 5);

  }

  int  s = 0;
  for (int i = 0; i < d->numBlocks(); i++) {
    REQUIRE(d->blockSize(i) == blockSizes[i]);
    REQUIRE(d->blockStart(i) == s);
    s += d->blockSize(i);
  }

  REQUIRE(d->numElements() == s);

  if (mat::DimensionDescriptor<B, NB>::Traits::numBlocksAtCompileTime == mat::Dynamic) {
    int b = 8;
    blockSizes.push_back(b);
    d->addBlock(b);


    s = 0;
    for (int i = 0; i < d->numBlocks(); i++) {
      REQUIRE(d->blockSize(i) == blockSizes[i]);
      REQUIRE(d->blockStart(i) == s);
      s += d->blockSize(i);
    }

    REQUIRE(d->numElements() == s);
  }
}

TEST_CASE("DimensionDescriptorTraits", "[DimensionDescriptor]") {

  {
    using T = mat::DimensionDescriptorTraits<3, 5>;
    static_assert(std::is_same<T::BlockSizeType, int>::value, "");
    static_assert(T::blockSizeAtCompileTime == 3, "");
    static_assert(T::blockType == mat::Fixed, "");
    static_assert(T::numBlocksAtCompileTime == 5, "");
    static_assert(T::numElementsAtCompileTime == 15, "");
  }

  {
    using T = mat::DimensionDescriptorTraits<3, mat::Dynamic>;
    static_assert(std::is_same<T::BlockSizeType, int>::value, "");
    static_assert(T::blockSizeAtCompileTime == 3, "");
    static_assert(T::blockType == mat::Fixed, "");
    static_assert(T::numBlocksAtCompileTime == mat::Dynamic, "");
    static_assert(T::numElementsAtCompileTime == mat::Dynamic, "");
  }

  {
    using T = mat::DimensionDescriptorTraits<mat::Dynamic, 5>;
    static_assert(std::is_same<T::BlockSizeType, int>::value, "");
    static_assert(T::blockSizeAtCompileTime == mat::Dynamic, "");
    static_assert(T::blockType == mat::Dynamic, "");
    static_assert(T::numBlocksAtCompileTime == 5, "");
    static_assert(T::numElementsAtCompileTime == mat::Dynamic, "");
  }

  {
    using T = mat::DimensionDescriptorTraits<mat::Dynamic, mat::Dynamic>;
    static_assert(std::is_same<T::BlockSizeType, int>::value, "");
    static_assert(T::blockSizeAtCompileTime == mat::Dynamic, "");
    static_assert(T::blockType == mat::Dynamic, "");
    static_assert(T::numBlocksAtCompileTime == mat::Dynamic, "");
    static_assert(T::numElementsAtCompileTime == mat::Dynamic, "");
  }


  {
    using T = mat::DimensionDescriptorTraits<mat::Variable, 5>;
    static_assert(std::is_same<T::BlockSizeType, std::vector<int>>::value, "this should be std::vector<int>");
    static_assert(T::blockSizeAtCompileTime == mat::Dynamic, "");
    static_assert(T::blockType == mat::Variable, "");
    static_assert(T::numBlocksAtCompileTime == 5, "");
    static_assert(T::numElementsAtCompileTime == mat::Dynamic, "");
  }

  {
    using T = mat::DimensionDescriptorTraits<mat::Variable, mat::Dynamic>;
    static_assert(std::is_same<T::BlockSizeType, std::vector<int>>::value, "this should be std::vector<int>");
    static_assert(T::blockSizeAtCompileTime == mat::Dynamic, "");
    static_assert(T::blockType == mat::Variable, "");
    static_assert(T::numBlocksAtCompileTime == mat::Dynamic, "");
    static_assert(T::numElementsAtCompileTime == mat::Dynamic, "");
  }


}
