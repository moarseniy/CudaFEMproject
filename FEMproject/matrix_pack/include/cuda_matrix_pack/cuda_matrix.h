#ifndef CUDA_MATRIX_H
#define CUDA_MATRIX_H

#include <iostream>
#include <cassert>

#include <matrix_pack/matrix_pack.h>
#include "cuda_matrix_pack/matrix_kernels.h"

class CUDA_Matrix : public Matrix {
private:
  static cublasHandle_t _handle;

  void _gpuClearData();
  void _gpuInitData();
public:
  static void _initCUDA();

  CUDA_Matrix();
  CUDA_Matrix(size_t numRows, size_t numCols);
  CUDA_Matrix(float *data, size_t numRows, size_t numCols, bool hasData = false);
  CUDA_Matrix(const CUDA_Matrix &like, bool copy = false); //TODO: make it const
  ~CUDA_Matrix();

//  void addWeighted();
  void product(Matrix &src, Matrix &tgt, bool a_tr = false, float scaleA = 1.f,
                                         bool b_tr = false, float scaleB = 1.f) override;
  void bmm(const Matrix &a, size_t subRowsA, size_t subColsA, bool trans_a,
           const Matrix &b, size_t subColsB, bool trans_b, size_t batchCount, const float alpha = 1.f) override;

  void copy(Matrix &target);
  void resize(size_t numRows, size_t numCols) override;
  void resize(const Matrix &like) override;
  void flatten() override;

  void divideElementwise(Matrix &v, Axis ax) override;
  void divideElementwise(Matrix &src, Matrix &target) override;
  void divideElementwise(Matrix &src) override;
  void divideElementwise(float value) override;
  void scale(float value) override;
  void scale(Matrix &v, Axis ax = ALL) override;
  float dotProduct(Matrix &src) override;
  void sort() override;
  void sort(Matrix &target) override;
  void sort_by_key(Matrix &keys) override;
  void reduce_by_key(Matrix &keys, Matrix &target) override;
  float det() override;
  float l2norm() override;
  void uniformRandomize(float v1 = 0.f, float v2 = 1.f) override;
  void fillSequence(float startValue) override;

  float min() override;
  float max() override;

  void setTo(float value) override;
  void getDiagonal(size_t size, Matrix &tgt) override;

//  std::unique_ptr<Matrix> subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol) const;
//  void subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol, Matrix &target) const;

  void multiplyByVec(const Matrix &vec, Matrix &target) const override;
  void addWeighted(Matrix &b, float alpha, float beta) override;
  void addWeighted(Matrix &b, float alpha, float beta, Matrix &target) override;
  void add(Matrix &src) override;
  void add(Matrix &src, Matrix &tgt) override;
  void subtract(Matrix &src) override;
  void subtract(Matrix &src, Matrix &tgt) override;
  static void copy(const Matrix &src, Matrix &tgt);
};


#endif // CUDA_MATRIX_H
