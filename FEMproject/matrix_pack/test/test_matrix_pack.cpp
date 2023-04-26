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

TEST(MathVectorOperations, reduce_by_key) {
  size_t rsz(2);
  size_t csz(4);
  {
    CUDA_Matrix values(rsz, csz), keys(rsz, csz), tgt(rsz, csz);

    values.uniformRandomize(0.f, 1.f);
    keys.fillSequence(0);

    // GPU reduction by key
    values.reduce_by_key(keys, tgt);

    CPU_Matrix hv, hk, hres(rsz, csz), res;
    values.copy(hv);
    keys.copy(hk);
    tgt.copy(res);

    // CPU reduction by key
    hv.reduce_by_key(hk, hres);

    float *gpu_data = res.get_data();
    float *cpu_data = hres.get_data();

    for (size_t i(0); i < res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
}

TEST(MatrixBinaryOperations, bmm) {
  size_t batchCount(10);
  size_t rsz(20);
  size_t csz(40);
  size_t msz(30);
  {
    CUDA_Matrix a(batchCount, csz * rsz), b(batchCount, csz * msz), c(batchCount, rsz * msz);

    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    // GPU Batched Matrix Multiplication (Row-major format)
    c.bmm(b, msz, csz, false, a, rsz, true, batchCount);

    CPU_Matrix ha, hb, hc, cpu_res(batchCount, rsz * msz);
    a.copy(ha);
    b.copy(hb);
    c.copy(hc);

    // CPU Batched Matrix Multiplication (Row-major format)
    cpu_res.bmm(hb, msz, csz, false, ha, rsz, true, batchCount);

    float *cdata = hc.get_data();
    float *cpu_data = cpu_res.get_data();

    for (size_t i(0); i < cpu_res.get_numElements(); ++i) {
      ASSERT_NEAR(cdata[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(batchCount, csz * rsz), b(batchCount, csz * msz), c(batchCount, rsz * msz);

    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    // GPU Batched Matrix Multiplication (Column-major format)
    c.bmm(a, rsz, csz, true, b, msz, false, batchCount);

    CPU_Matrix ha, hb, hc, cpu_res(batchCount, rsz * msz);
    a.copy(ha);
    b.copy(hb);
    c.copy(hc);

    // CPU Batched Matrix Multiplication (Column-major format)
    cpu_res.bmm(ha, rsz, csz, true, hb, msz, false, batchCount);

    float *cdata = hc.get_data();
    float *cpu_data = cpu_res.get_data();

    for (size_t i(0); i < cpu_res.get_numElements(); ++i) {
      ASSERT_NEAR(cdata[i], cpu_data[i], EPS);
    }
  }
}

TEST(MatrixBinaryOperations, addWeighted) {
  size_t rsz(200);
  size_t csz(400);
  {
    CUDA_Matrix a(rsz, csz), b(rsz, csz), c(rsz, csz);

    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    // GPU Weighted Addition
    a.addWeighted(b, 0.5f, 0.1f, c);

    CPU_Matrix ha, hb, hc, cpu_res;
    a.copy(ha);
    b.copy(hb);
    c.copy(hc);

    // CPU Weighted Addition
    ha.addWeighted(hb, 0.5f, 0.1f, cpu_res);

    float *cdata = hc.get_data();
    float *cpu_data = cpu_res.get_data();

    for (size_t i(0); i < cpu_res.get_numElements(); ++i) {
      ASSERT_NEAR(cdata[i], cpu_data[i], EPS);
    }
  }
}

TEST(MathVectorOperations, l2norm_diff) {
  size_t sz(200);
  {
    CUDA_Matrix a(1, sz), b(1, sz), c(1, sz);

    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    // GPU L2 Euclidean norm of difference
    a.subtract(b, c);
    float gpu_res = c.l2norm();

    CPU_Matrix ha, hb, hc;
    a.copy(ha);
    b.copy(hb);

    // CPU L2 Euclidean norm of difference
    ha.subtract(hb, hc);
    float cpu_res = c.l2norm();

    ASSERT_NEAR(cpu_res, gpu_res, EPS);
  }
}

TEST(MatrixOperations, sort) {
  size_t rsz(4);
  size_t csz(3);
  {
    CUDA_Matrix a(rsz, csz), b;
    a.uniformRandomize(0.f, 1.f);

    a.sort(b);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(res);

    ha.sort(hb);

    float *gpu_data = res.get_data();
    float *cpu_data = hb.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_EQ(gpu_data[i], cpu_data[i]);
    }
  }
}

TEST(MatrixBinaryOperations, scale) {
  size_t rsz(40);
  size_t csz(30);
  {
    CUDA_Matrix a(rsz, csz), b(1, rsz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.scale(b, X);
    a.copy(res);

    ha.scale(hb, X);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz), b(1, csz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.scale(b, Y);
    a.copy(res);

    ha.scale(hb, Y);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz), b(rsz, csz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.scale(b, ALL);
    a.copy(res);

    ha.scale(hb, ALL);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz);
    a.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, res;
    a.copy(ha);

    a.scale(0.66f);
    a.copy(res);

    ha.scale(0.66f);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < ha.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
}

TEST(MatrixBinaryOperations, divide) {
  size_t rsz(40);
  size_t csz(30);
  {
    CUDA_Matrix a(rsz, csz), b(1, rsz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.divideElementwise(b, X);
    a.copy(res);

    ha.divideElementwise(hb, X);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz), b(1, csz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.divideElementwise(b, Y);
    a.copy(res);

    ha.divideElementwise(hb, Y);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz), b(rsz, csz);
    a.uniformRandomize(0.f, 1.f);
    b.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, hb, res;
    a.copy(ha);
    b.copy(hb);

    a.divideElementwise(b, ALL);
    a.copy(res);

    ha.divideElementwise(hb, ALL);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < hb.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
  {
    CUDA_Matrix a(rsz, csz);
    a.uniformRandomize(0.f, 1.f);

    CPU_Matrix ha, res;
    a.copy(ha);

    a.divideElementwise(0.66f);
    a.copy(res);

    ha.divideElementwise(0.66f);

    float *gpu_data = res.get_data();
    float *cpu_data = ha.get_data();

    for (size_t i(0); i < ha.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS);
    }
  }
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

//  std::cout << "\nRes\n";
//  for (auto &key : keys) {
//    auto range = mymap.equal_range(key);
//    std::cout << std::accumulate(range.first, range.second,
//        0.0f,
//        [](int a, std::pair<int, int> b) { return a + b.second; }) << " ";
//  }

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

int main(int argc, char *argv[]) {
  std::srand(time(0));
  CUDA_Matrix::_initCUDA();

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
