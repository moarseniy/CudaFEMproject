#ifndef MATRIX_PACK_H
#define MATRIX_PACK_H

#include <cassert>
#include <stdexcept>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>

#include <Tools.h>

enum DEVICE_NAME {
  CPU,
  CUDA,
  ONEAPI
};

enum Axis {
  X,
  Y,
  ALL
};

class Matrix {
public:
  Matrix(size_t numRows, size_t numCols, bool hasData, float *data, DEVICE_NAME device);
  virtual ~Matrix() {};

  void setDevice(DEVICE_NAME device);

//  virtual void addWeighted() = 0;
  virtual void product(Matrix &src, Matrix &tgt, bool a_tr = false, float scaleA = 1.f,
                                                 bool b_tr = false, float scaleB = 1.f) = 0;
  virtual void bmm(const Matrix &a, size_t subRowsA, size_t subColsA, bool trans_a,
                   const Matrix &b, size_t subColsB, bool trans_b, size_t batchCount, const float alpha = 1.f) = 0;

  virtual void divideElementwise(Matrix &v, Axis ax) = 0;
  virtual void divideElementwise(Matrix &src, Matrix &tgt) = 0;
  virtual void divideElementwise(Matrix &src) = 0;
  virtual void divideElementwise(float value) = 0;
  virtual void scale(float value) = 0;
  virtual void scale(Matrix &v, Axis ax = ALL) = 0;
  virtual float dotProduct(Matrix &src) = 0;
  virtual void sort(Matrix &target) = 0;
  virtual void sort() = 0;
  virtual void sort_by_key(Matrix &keys) = 0;
  virtual void reduce_by_key(Matrix &keys, Matrix &target) = 0;
  virtual float det() = 0;
  virtual float l2norm() = 0;

  // Right-vector multiplication
  virtual void multiplyByVec(const Matrix &vec, Matrix &target) const = 0;

  virtual void addWeighted(Matrix &b, float alpha, float beta) = 0;
  virtual void addWeighted(Matrix &b, float alpha, float beta, Matrix &target) = 0;

  virtual void add(Matrix &src) = 0;
  virtual void add(Matrix &src, Matrix &tgt) = 0;
  virtual void subtract(Matrix &src) = 0;
  virtual void subtract(Matrix &src, Matrix &tgt) = 0;

  void copy(Matrix &tgt) const;
  virtual void resize(size_t numRows, size_t numCols) = 0;
  virtual void resize(const Matrix &like) = 0;
  virtual void flatten() = 0;

  static Matrix* setMatrix(DEVICE_NAME device, size_t numRows, size_t numCols);
  static Matrix* setMatrix(DEVICE_NAME device);
  static Matrix* setMatrix(Matrix &like);

  static Matrix* setVector(DEVICE_NAME device, size_t size);
  static Matrix* setVector(DEVICE_NAME device);
  static Matrix* setVector(Matrix &like);

  virtual void setTo(float value) = 0;

  std::unique_ptr<Matrix> subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol) const;
  void subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol, Matrix &target) const;

  std::unique_ptr<Matrix> getRows(size_t startRow, size_t endRow) const;
  void getRows(size_t startRow, size_t endRow, Matrix &target) const;

  std::unique_ptr<Matrix> getRow(size_t numRow) const;
  void getRow(size_t numRow, Matrix &target) const;

  void Show();
  void writeToFile(const std::string path);
  bool isSameAs(Matrix &src) const;

  DEVICE_NAME get_device() const;
  size_t get_numRows() const;
  size_t get_numCols() const;
  size_t get_numElements() const;
  float *get_data() const;

  float& operator [](size_t index);
  float& operator ()(size_t i, size_t j);

  virtual void uniformRandomize(float v1 = 0.f, float v2 = 1.f) = 0;
  virtual void fillSequence(float startValue) = 0;

protected:

  DEVICE_NAME _device;

  size_t _numRows;
  size_t _numCols;
  size_t _numElements;
  float *_data;

  bool _hasData;
};

class CPU_Matrix : public Matrix {
public:
  CPU_Matrix();
  CPU_Matrix(size_t numRows, size_t numCols);
  CPU_Matrix(float *data, size_t numRows, size_t numCols, bool hasData = false);
  CPU_Matrix(const CPU_Matrix &like, bool copy = false);
  ~CPU_Matrix();

  // A =
//  void addWeighted();
  void product(Matrix &src, Matrix &tgt, bool a_tr = false, float scaleA = 1.f,
                                         bool b_tr = false, float scaleB = 1.f) override;

  void bmm(const Matrix &a, size_t subRowsA, size_t subColsA, bool trans_a,
           const Matrix &b, size_t subColsB, bool trans_b, size_t batchCount, const float alpha = 1.f) override;

  void addWeighted(Matrix &b, float alpha, float beta) override;
  void addWeighted(Matrix &b, float alpha, float beta, Matrix &target) override;
  void subtract(Matrix &src) override;
  void subtract(Matrix &src, Matrix &tgt) override;
  void add(Matrix &src) override;
  void add(Matrix &src, Matrix &tgt) override;

  void copy(Matrix &tgt);
  void resize(size_t numRows, size_t numCols) override;
  void resize(const Matrix &like) override;
  void flatten() override;

  static void copy(const Matrix &src, Matrix &tgt);

  void divideElementwise(Matrix &v, Axis) override;
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

  void setTo(float value) override;
//  std::unique_ptr<Matrix> subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol) const;
//  void subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol, Matrix &target) const;

  void multiplyByVec(const Matrix &vec, Matrix &target) const override;

//  float& operator [](size_t index);
//  float& operator ()(size_t i, size_t j);
};

#endif // MATRIX_PACK_H
