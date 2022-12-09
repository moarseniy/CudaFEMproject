#include <matrix_pack/matrix_pack.h>
#include <matrix_pack_cuda/gpu_matrix.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/compute/system.hpp>
#include <boost/compute/algorithm/fill.hpp>
#include <boost/compute/algorithm/reduce_by_key.hpp>
#include <boost/compute/container/vector.hpp>

TEST(MatrixFuncs, smth) {
//  const int N = 7;
//  int *A = new int[N] {1, 3, 3, 3, 2, 2, 1}; // input keys
//  int B[N] = {9, 8, 7, 6, 5, 4, 3}; // input values
//  int C[N];                         // output keys
//  int D[N];                         // output values
//  boost::compute::reduce_by_key(A, A + N, B, C, D);

  EXPECT_EQ(true, true);
}

int main(int argc, char *argv[]) {
  CPU_Matrix A(7, 10);

//  GPU_Matrix c(5, 2);
  Matrix *c = new CUDA_Matrix(7, 10);
  CUDA_Matrix d(7, 10);

//  c->copy(d, true);
  c->copy(A, false);
  A.Show();
  delete c;
//  delete B;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
