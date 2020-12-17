#include <Eigen/Core>
#include <vector>
#include <chrono>
#include <iostream>
#include <Eigen/Sparse>

void testBlocks() {
  // observe in debug what a block contains
  Eigen::Matrix<double, 10, 20> ms;
  Eigen::MatrixXd md(10, 20);

  auto bs = ms.block<2, 2>(0, 0, 2, 2);
  auto bm = md.block<2, 2>(0, 0, 2, 2);

  auto bbs = bs.block<2, 1>(0, 0, 2, 1);

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

void testMapWithStride() {
  {
    Eigen::VectorXf v(18);
    for (int i = 0; i < v.size(); i++) {
      v(i) = i;
    }

    Eigen::Map<Eigen::Matrix3f, 0, Eigen::OuterStride<>> m33(v.data(), 3, 3, Eigen::OuterStride<>(6));
    std::cout << m33 << std::endl << std::endl;

  }
  // sparse mat
  {
    Eigen::SparseMatrix<float> s(8,4);
    s.reserve(20);

    int c = 0;

    //column 0
    s.startVec(0);
    s.insertBack(0, 0) = c++;
    s.insertBack(1, 0) = c++;

    s.insertBack(4, 0) = c++;
    s.insertBack(5, 0) = c++;

    //column 1
    s.startVec(1);
    s.insertBack(0, 1) = c++;
    s.insertBack(1, 1) = c++;

    s.insertBack(4, 1) = c++;
    s.insertBack(5, 1) = c++;

    //column 2
    s.startVec(2);
    s.insertBack(0, 2) = c++;
    s.insertBack(1, 2) = c++;

    s.insertBack(2, 2) = c++;
    s.insertBack(3, 2) = c++;

    s.insertBack(4, 2) = c++;
    s.insertBack(5, 2) = c++;

    //column 3
    s.startVec(3);
    s.insertBack(0, 3) = c++;
    s.insertBack(1, 3) = c++;

    s.insertBack(2, 3) = c++;
    s.insertBack(3, 3) = c++;

    s.insertBack(4, 3) = c++;
    s.insertBack(5, 3) = c++;

    s.finalize();

    Eigen::MatrixXf d = s;


    //Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> mXX(nullptr);

    Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> m00(s.coeffs().data(), 2, 2, Eigen::OuterStride<>(4));
    std::cout << "B_00 " << std::endl << d.block<2, 2>(0, 0) << std::endl << std::endl;
    std::cout << "m_00 " << std::endl << m00 << std::endl << std::endl;
    
    Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> m20(s.coeffs().data() + 2, 2, 2, Eigen::OuterStride<>(4));
    std::cout << "B_20 " << std::endl << d.block<2, 2>(4, 0) << std::endl << std::endl;
    std::cout << "m_20 " << std::endl << m20 << std::endl << std::endl;

    Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> m01(s.coeffs().data() + 8, 2, 2, Eigen::OuterStride<>(6));
    std::cout << "B_01 " << std::endl << d.block<2, 2>(0, 2) << std::endl << std::endl;
    std::cout << "m_01 " << std::endl << m01 << std::endl << std::endl;

    Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> m11(s.coeffs().data() + 10, 2, 2, Eigen::OuterStride<>(6));
    std::cout << "B_11 " << std::endl << d.block<2, 2>(2, 2) << std::endl << std::endl;
    std::cout << "m_11 " << std::endl << m11 << std::endl << std::endl;

    Eigen::Map<Eigen::Matrix2f, 0, Eigen::OuterStride<>> m21(s.coeffs().data() + 12, 2, 2, Eigen::OuterStride<>(6));
    std::cout << "B_21 " << std::endl << d.block<2, 2>(4, 2) << std::endl << std::endl;
    std::cout << "m_21 " << std::endl << m21 << std::endl << std::endl;

    m21 = m21 * 2;
    std::cout << "B_21 " << std::endl << d.block<2, 2>(4, 2) << std::endl << std::endl;
    std::cout << "m_21 " << std::endl << m21 << std::endl << std::endl;

  }

}

int main(int argc, char* argv[]) {
  // what is inside Block?
  testBlocks();

  testMapWithStride();


  // slicing

  // timing 
  coutExecTime(testSetZeroTime);
  coutExecTime(testSetZeroTimeParallel);
  coutExecTime(testSetZeroSingleMatTime);
  coutExecTime(testSetZeroTimeDyn1);
  coutExecTime(testSetZeroTimeDyn2);

  return 0;
}