#ifndef LINAL_H
#define LINAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <stdio.h>
#include <cmath>

#include "Tools.h"

class MyArray;
class Matrix;

class MyArray {
public:
  MyArray();
  MyArray(int array_size);
  MyArray(const MyArray &a);
  MyArray(int array_size, float a, float b);
  ~MyArray();
  void Resize(int new_size);
  void Show();
  void ShowNonzero();
  void zap();
  void Set(int index, float value);
  void add(MyArray a);
  void add_weighted(MyArray a, float v1, float v2);
  void scale(float value);
  int get_size();
  float* get_data();
  float & operator [](int index);
  MyArray operator =(const MyArray &a);
  float norma();
  float CNorm();
  float dot_product(MyArray v2);
  bool equalsToArray(MyArray a, float eps);
  MyArray multiplyElementwise(MyArray a);
  MyArray divideByElementwise(MyArray a);
  //    void multiplyElementwise(MyArray a);
  //    void divideByElementwise(MyArray a);
  //    void multiplyElementwise(MyArray a, MyArray &res);
  //    void divideByElementwise(MyArray a, MyArray &res);

  void grad(float (*f)(MyArray x), MyArray &res, float eps);
  void Hessian(float (*f)(MyArray x), Matrix &matr, float eps);

  void WriteToFile();
  void ReadFromFile();
private:
  int array_size;
  float *p;
};

class Matrix{
public:
  Matrix();
  Matrix(int n);
  Matrix(int row, int col);
  Matrix(const Matrix &a);
  Matrix(int row, int col, bool isDiag);
  ~Matrix();
  void Resize(int row, int col);
  void Show();
  float & operator ()(int i,int j);
  void Set(int index1,int index2,float value);
  float Get(int index1,int index2);
  int get_row();
  int get_col();
  float* get_data();
  void zap();
  Matrix operator =(const Matrix &a);
  void LU_decomposition(Matrix &L, Matrix &U, int n);
  void Solve_Gauss_reverse(Matrix matr, MyArray B, MyArray &res, int n, bool isDown);
  void LU_solve(MyArray B, MyArray &result, int n);
  void CGM_solve(Matrix A, MyArray B, MyArray &result, int n);
  void scale(float value);
  Matrix Sum(Matrix &a);
  Matrix weightedSum(Matrix &a, float alpha, float beta);
  Matrix Difference(Matrix &a);
  Matrix Product(Matrix &a);
  MyArray Product(MyArray &v);
  int CountNonzero();
  Matrix & transpose();
  bool equalsToMatrix(Matrix &a, float eps);

  void inverse3x3(Matrix &res);

  void transpose2();
  void transpose_not_square();

  Matrix & Gauss();
  float det_gauss();
  float det(int n);
  void Get_matrix(int n, Matrix &temp_matr, int indRow, int indCol);
  void inverse(Matrix &matr, int n, bool isDiag);

  void WriteToFile();
  void ReadFromFile();

private:
  int row;
  int col;
  bool isDiag;
  bool isSym;
  float *m;
};

#endif
