#pragma once

#include <Eigen/Core>
#include <vector>
#include <chrono>
#include <iostream>

void testBlocks() {
  // observe in debug what a block contains
  Eigen::Matrix<double, 10, 20> ms;
  Eigen::MatrixXd md(10, 20);

  auto bs = ms.block<2, 2>(0, 0, 2, 2);
  auto bm = md.block<2, 2>(0, 0, 2, 2);

}

const int N = 1000000;
const int S = 3;

void testSetZeroTime() {
  std::cout << __FUNCTION__ << std::endl;

  auto t0 = std::chrono::steady_clock::now();
  std::vector<Eigen::Matrix<double, S, S>> m(N);
  auto t1 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m[i].setZero();
  }
  auto t2 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t0).count() << std::endl;

}

void testSetZeroTimeDyn() {
  std::cout << __FUNCTION__ << std::endl;
  {
    auto t0 = std::chrono::steady_clock::now();
    std::vector<Eigen::MatrixXd> m(N, Eigen::MatrixXd(S, S));
    auto t1 = std::chrono::steady_clock::now();

    for (int i = 0; i < N; i++) {
      m[i].setZero();
    }
    auto t2 = std::chrono::steady_clock::now();

    std::cout << "creation (copy)" << std::chrono::duration<double>(t1 - t0).count() << std::endl;
    std::cout << "zero " << std::chrono::duration<double>(t2 - t0).count() << std::endl;
  }
  {
    auto t0 = std::chrono::steady_clock::now();
    std::vector<Eigen::MatrixXd> m(N);
    for (int i = 0; i < N; i++) {
      m[i].resize(S, S);
    }
    auto t1 = std::chrono::steady_clock::now();

    for (int i = 0; i < N; i++) {
      m[i].setZero();
    }
    auto t2 = std::chrono::steady_clock::now();

    std::cout << "creation (loop)" << std::chrono::duration<double>(t1 - t0).count() << std::endl;
    std::cout << "zero " << std::chrono::duration<double>(t2 - t0).count() << std::endl;
  }
}

void testSetZeroSingleMatTime() {
  std::cout << __FUNCTION__ << std::endl;

  auto t0 = std::chrono::steady_clock::now();
  Eigen::Matrix<double, S, Eigen::Dynamic> m(S, N*S);
  auto t1 = std::chrono::steady_clock::now();

  m.setZero();
  
  auto t2 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t0).count() << std::endl;

}

int main(int argc, char* argv[]) {

  testBlocks();

  testSetZeroTime();
  testSetZeroTimeDyn();
  testSetZeroSingleMatTime();

  return 0;
}