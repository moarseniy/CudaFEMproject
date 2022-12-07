#include <matrix_pack/matrix_pack.h>
#include <matrix_pack_cuda/gpu_matrix.h>

//#include <gtest/gtest.h>

//TEST(MatrixFuncs, smth) {

//}

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
//  ::testing::InitGoogleTest(&argc, argv);
//  return RUN_ALL_TESTS();
}
