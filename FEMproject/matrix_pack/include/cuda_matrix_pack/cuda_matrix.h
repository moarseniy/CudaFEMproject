#ifndef CUDA_MATRIX_H
#define CUDA_MATRIX_H

#include <iostream>
#include <cassert>

//#include <cublas.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

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
#include "cuda_matrix_pack/matrix_kernels.h"

class Matrix;
class CUDA_Matrix : public Matrix {
private:
  void gpuClearData();
  void gpuInitData();
public:
  CUDA_Matrix();
  CUDA_Matrix(size_t numRows, size_t numCols);
  CUDA_Matrix(float *data, size_t numRows, size_t numCols, bool hasData = false);
  CUDA_Matrix(const CUDA_Matrix &like);
  ~CUDA_Matrix();

//  void addWeighted();
  void product(Matrix &src, Matrix &tgt) override;

  void add(Matrix &src) override;

  void copy(Matrix &target);
  void resize(size_t numRows, size_t numCols) override;
  void divideElementwise(Matrix &src) override;
  void divideElementwise(float value) override;
  void scale(float value) override;
  float dotProduct(Matrix &src) override;
  void sort() override;
  void sort_by_key(Matrix &keys) override;
  void reduce_by_key(Matrix &keys, Matrix &target) override;
  float det() override;

  void setTo(float value) override;
//  std::unique_ptr<Matrix> subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol) const;
//  void subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol, Matrix &target) const;

  void multiplyByVec(const Matrix &vec, Matrix &target) const override;
  void addWeighted(Matrix &b, float alpha, float beta) override;
  static void copy(Matrix &src, Matrix &tgt);
};


#endif // CUDA_MATRIX_H
