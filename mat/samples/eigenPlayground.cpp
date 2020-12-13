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

std::chrono::steady_clock::time_point testSetZeroTime() {
  std::cout << std::endl << "----" << __FUNCTION__ << std::endl;

  auto t0 = std::chrono::steady_clock::now();
  std::vector<Eigen::Matrix<double, S, S>> m(N);
  auto t1 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m[i].setZero();
  }
  auto t2 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m[i].setIdentity();
  }
  auto t3 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
  std::cout << "setIdentity x block " << std::chrono::duration<double>(t3 - t2).count() << std::endl;
  std::cout << "TOT " << std::chrono::duration<double>(t3 - t0).count() << std::endl;

  return t3;
}

std::chrono::steady_clock::time_point testSetZeroTimeParallel() {
  std::cout << std::endl << "----" << __FUNCTION__ << std::endl;

  auto t0 = std::chrono::steady_clock::now();
  std::vector<Eigen::Matrix<double, S, S>> m(N);
  auto t1 = std::chrono::steady_clock::now();

#pragma omp parallel for default (shared) schedule(static)
  for (int i = 0; i < N; i++) {
    m[i].setZero();
  }
  auto t2 = std::chrono::steady_clock::now();

#pragma omp parallel for default (shared) schedule(static)
  for (int i = 0; i < N; i++) {
    m[i].setIdentity();
  }
  auto t3 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
  std::cout << "setIdentity x block " << std::chrono::duration<double>(t3 - t2).count() << std::endl;
  std::cout << "TOT " << std::chrono::duration<double>(t3 - t0).count() << std::endl;
  return t3;
}


std::chrono::steady_clock::time_point testSetZeroTimeDyn1() {
  std::cout << std::endl << "----" << __FUNCTION__ << std::endl;
  auto t0 = std::chrono::steady_clock::now();
  std::vector<Eigen::MatrixXd> m(N, Eigen::MatrixXd(S, S));
  auto t1 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m[i].setZero();
  }
  auto t2 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m[i].setIdentity();
  }
  auto t3 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
  std::cout << "setIdentity x block " << std::chrono::duration<double>(t3 - t2).count() << std::endl;
  std::cout << "TOT " << std::chrono::duration<double>(t3 - t0).count() << std::endl;
  return t3;
}

std::chrono::steady_clock::time_point testSetZeroTimeDyn2() {
  std::cout << std::endl << "----" << __FUNCTION__ << std::endl;
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

  for (int i = 0; i < N; i++) {
    m[i].setIdentity();
  }
  auto t3 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
  std::cout << "setIdentity x block " << std::chrono::duration<double>(t3 - t2).count() << std::endl;
  std::cout << "TOT " << std::chrono::duration<double>(t3 - t0).count() << std::endl;
  return t3;
}

std::chrono::steady_clock::time_point testSetZeroSingleMatTime() {
  std::cout << std::endl << "----" << __FUNCTION__ << std::endl;

  auto t0 = std::chrono::steady_clock::now();
  Eigen::Matrix<double, S, Eigen::Dynamic> m(S, N * S);
  auto t1 = std::chrono::steady_clock::now();

  m.setZero();

  auto t2 = std::chrono::steady_clock::now();

  for (int i = 0; i < N; i++) {
    m.block<S, S>(0, i * S, S, S).setIdentity();
  }
  auto t3 = std::chrono::steady_clock::now();

#pragma omp parallel for default (shared) schedule(static)
  for (int i = 0; i < N; i++) {
    m.block<S, S>(0, i * S, S, S).setIdentity();
  }
  auto t4 = std::chrono::steady_clock::now();

  std::cout << "creation " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  std::cout << "zero " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
  std::cout << "setIdentity x block " << std::chrono::duration<double>(t3 - t2).count() << std::endl;
  std::cout << "TOT " << std::chrono::duration<double>(t3 - t0).count() << std::endl;

  std::cout << "setIdentity x block (parallel)" << std::chrono::duration<double>(t4 - t3).count() << std::endl;
  return t3;
}

void coutExecTime(std::chrono::steady_clock::time_point(*f)()) {
  auto t0 = f();
  auto t1 = std::chrono::steady_clock::now();
  std::cout << "-- exit time " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
}

int main(int argc, char* argv[]) {

  testBlocks();

  // timing 
  coutExecTime(testSetZeroTime);
  coutExecTime(testSetZeroTimeParallel);
  coutExecTime(testSetZeroSingleMatTime);
  coutExecTime(testSetZeroTimeDyn1);
  coutExecTime(testSetZeroTimeDyn2);

  return 0;
}