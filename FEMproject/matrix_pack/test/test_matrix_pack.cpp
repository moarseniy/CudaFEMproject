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

#define EPS 1e-4

TEST(MatrixFuncs, sort) {
  CPU_Matrix m(4, 4);
  m.fillRandValues(0.f, 1.f);
  m.Show();
  CUDA_Matrix m_d;
  m.copy(m_d);

  m_d.sort();
  exit(-1);
}

TEST(MatrixFuncs, Sdgmm) {
  CPU_Matrix A(2, 4), v(4, 1);
  CUDA_Matrix A_d, v_d;

  A.fillRandValues(0.f, 5.f);
  A.copy(A_d);
  v.fillRandValues(0.f, 5.f);
  v.copy(v_d);
  A.Show();
   std::cout << "\n";
  v.Show();
   std::cout << "\n";

  A.scale(v, Y);
  A.Show();
  std::cout << "\n";

  A_d.scale(v_d, Y);
  CPU_Matrix tmp;
  A_d.copy(tmp);
  tmp.Show();
  exit(-1);
}

TEST(MatrixFuncs, Sgemm) {
  Matrix *A = Matrix::setMatrix(CPU, 1, 3 * 6);
  for (size_t i = 0; i < A->get_numElements(); ++i)
    (*A)[i] =i;
  std::cout << "A\n";
  A->Show();
  Matrix *B = Matrix::setMatrix(CPU, 1, 3 * 3);
  for (size_t i = 0; i < B->get_numElements(); ++i)
    (*B)[i] =i;
  std::cout << "B\n";
  B->Show();

  Matrix *C = Matrix::setMatrix(CPU, 1, 3 * 6);
  A->multiplyByVec(*B, *C);
  std::cout << "C\n";
  C->Show();

  Matrix *A_d = Matrix::setMatrix(CUDA, A->get_numRows(), A->get_numCols());
  Matrix *B_d = Matrix::setMatrix(CUDA, B->get_numRows(), B->get_numCols());
  Matrix *C_d = Matrix::setMatrix(CUDA, C->get_numRows(), C->get_numCols());

  C_d->setTo(0.f);
  A->copy(*A_d);
  B->copy(*B_d);



  C->setTo(0.f);
//  B->flatten();
//  C->bmm(*B, 2, 2, false, *A, 3, true, 1); // (A(2x3))^T * B(2x2) = C(3x2)
  C->bmm(*B, 3, 3, false, *A, 6, true, 1); // A(3x2) * B(2x3) = C(3x3)
  std::cout << "C\n";
  C->Show();

  C_d->bmm(*B_d, 3, 3, false, *A_d, 6, true, 1);
  CPU_Matrix tmp;
  C_d->copy(tmp);
//////  tmp.copy(*C_d);
  std::cout << "C\n";
  tmp.Show();
  exit(-1);
}

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

TEST(MatrixFuncs, mult) {
  std::cout << "product\n";

  CPU_Matrix m1(2, 3);
  m1.fillRandValues(0.f, 1.f);
  m1.Show();

  CPU_Matrix m2(2, 2);
  m2.fillRandValues(0.f, 1.f);
  m2.Show();

  CPU_Matrix m3(3, 2);
  m1.product(m2, m3, true);
  m3.Show();
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
