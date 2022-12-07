#ifndef CUDA_MATRIX_H
#define CUDA_MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <cuda_runtime.h>
#include <cuda.h>

#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/inner_product.h>

#include <matrix_pack/matrix_pack.h>
#include "matrix_pack_cuda/matrix_kernels.h"

class CUDA_Matrix : public Matrix {
private:
  void gpuClearData();
  void gpuInitData();
public:
  CUDA_Matrix();
  CUDA_Matrix(size_t numRows, size_t numCols);
  CUDA_Matrix(const CUDA_Matrix &like);
  ~CUDA_Matrix();

//  void addWeighted();
  void product(Matrix &src, Matrix &tgt) override;

  void copy(Matrix &target, bool isGPU) const override;
  void divideElementwise(Matrix &src) override;
  float dotProduct(Matrix &src) override;
  void sort() override;
  void sort_by_key(Matrix &keys) override;
  void reduce_by_key(Matrix &keys, Matrix &target) override;
  void multiplyByVec(const Matrix &vec, Matrix &target) const override;
};


#endif // CUDA_MATRIX_H
