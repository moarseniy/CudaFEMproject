#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

#include <gtest/gtest.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include <unordered_set>
#include <map>


TEST(MatrixFuncs, reduce) {
  const int N = 7;
  int A[N] = {1, 1, 3, 3, 3, 2, 2}; // input keys
  int B[N] = {9, 3, 8, 7, 6, 5, 4}; // input values
  int C[N];                         // output keys
  int D[N];                         // output values
//  thrust::reduce_by_key(thrust::host, A, A + N, B, C, D);
//  for (int i = 0; i < N; ++i) {
//    std::cout << D[i] << " ";
//  }
//  for (int i = 0; i < N; ++i) {
//    std::cout << C[i] << " ";
//  }
}

TEST(MatrixFuncs, smth) {
  const int N = 7;

  int *A = new int[N] {1, 1, 3, 3, 3, 2, 2}; // input keys
  int *B = new int[N] {9, 3, 8, 7, 6, 5, 4}; // input values
  int *C = new int[N];                         // output keys
  int *D = new int[N];                         // output values

  std::multimap<int, int> mymap;
  std::set<int> keys;
  for (int i = 0; i < N; ++i) {
    mymap.insert(std::multimap<int, int>::value_type(A[i], B[i]));
    keys.insert(A[i]);
  }

  std::cout << "\nRes\n";
  for (auto &key : keys) {
    auto range = mymap.equal_range(key);
    std::cout << std::accumulate(range.first, range.second,
        0.0f,
        [](int a, std::pair<int, int> b) { return a + b.second; }) << " ";
  }

  CPU_Matrix test1(1, 5), k(1, 5);
  test1.get_data()[0] = 1;
  test1.get_data()[1] = 2;
  test1.get_data()[2] = 3;
  test1.get_data()[3] = 4;
  test1.get_data()[4] = 5;
  k.get_data()[0] = 5;
  k.get_data()[1] = 5.5;
  k.get_data()[2] = 3;
  k.get_data()[3] = -1;
  k.get_data()[4] = 3;
  test1.sort_by_key(k);
//  test1.Show();
//  k.Show();

  EXPECT_EQ(true, true);
}

TEST(MatrixFuncs, copy) {
  CPU_Matrix m1(4, 4);

  m1.fillRandValues(0.f, 1.f);
  m1.Show();
  std::cout << "getRows\n";
  std::unique_ptr<Matrix> m3 = m1.getRows(1, 3);
  m3->setTo(1.f);
  m1.Show();
}

int main(int argc, char *argv[]) {
  CPU_Matrix A(7, 10);

//  GPU_Matrix c(5, 2);
  Matrix *c = new CUDA_Matrix(7, 10);
  CUDA_Matrix d(7, 10);

//  A.Show();
  delete c;
//  delete B;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
