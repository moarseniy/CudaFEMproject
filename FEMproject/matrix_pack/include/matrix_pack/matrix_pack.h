#ifndef MATRIX_PACK_H
#define MATRIX_PACK_H

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <cublas.h>
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

class Matrix {
public:
  Matrix(size_t numRows, size_t numCols, float *data);
  virtual ~Matrix() {};

//  virtual void addWeighted() = 0;
  virtual void product(Matrix &src, Matrix &tgt) = 0;

  virtual void divideElementwise(Matrix &src) = 0;
  virtual float dotProduct(Matrix &src) = 0;
  virtual void sort() = 0;
  virtual void sort_by_key(Matrix &keys) = 0;
  virtual void reduce_by_key(Matrix &keys, Matrix &target) = 0;

  // Right-vector multiplication
  virtual void multiplyByVec(const Matrix &vec, Matrix &target) const = 0;

  virtual void copy(Matrix &target, bool isGPU) const = 0;
  void Show();

  size_t get_numRows() const;
  size_t get_numCols() const;
  size_t get_numElements() const;
  float *get_data() const;

  void fillRandValues(float v1, float v2);

protected:
  size_t _numRows;
  size_t _numCols;
  size_t _numElements;
  float *_data;
};

class CPU_Matrix : public Matrix {
public:
  CPU_Matrix();
  CPU_Matrix(size_t numRows, size_t numCols);
  CPU_Matrix(const CPU_Matrix &like);
  ~CPU_Matrix();

  // A =
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

#endif // MATRIX_PACK_H
