#ifndef LINAL_H
#define LINAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

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
    float dot_product(MyArray v2);
    bool equalsToArray(MyArray a, float eps);
    MyArray multiplyElementwise(MyArray a);
    MyArray divideByElementwise(MyArray a);

	void grad(float (*f)(MyArray x), MyArray &res, float eps);
	void Hessian(float (*f)(MyArray x), Matrix &matr, float eps);

	void WriteToFile();
    void ReadFromFile();
private:
	int array_size;
	float *p;
};

struct Triplet {
    Triplet(int x_value, int y_value, float value) {
        this->x_value = x_value;
        this->y_value = y_value;
        this->value = value;
    }
    int get_x() { return x_value; }
    int get_y() { return y_value; }
    float get_value() { return value; }

    void Show() {
        std::cout << x_value << " " << y_value << " " << value << "\n";
    }
    int x_value;
    int y_value;
    float value;
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
	Matrix Difference(Matrix &a);
	Matrix Product(Matrix &a);
    MyArray Product(MyArray &v);
	Matrix & transpose();
    bool equalsToMatrix(Matrix &a, float eps);

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

struct couple{
    couple(int x, int y){
        this->x = x;
        this->y = y;
    }
    int x;
    int y;
};

class SparseMatrixCOO {
public:
    SparseMatrixCOO();
    SparseMatrixCOO(int sparse_size);
    SparseMatrixCOO(const SparseMatrixCOO &m);
    SparseMatrixCOO operator =(const SparseMatrixCOO &m);
	~SparseMatrixCOO();
    void resize();
    void resize(int size);
    void ConvertTripletToSparse(std::vector<Triplet> t);
    void ConvertToMatrix(Matrix& M);
    void SortIt();
    void ConvertToCSR(int *ptr, int *ind, float *data_csr, int n);
    void SparseLU();
    void CGM_solve(MyArray B, MyArray &x_k, float eps);
    void PCG_solve(MyArray B, MyArray &x_k, float eps);
    MyArray MyltiplyByVector(MyArray v);
    MyArray MyltiplyTransposedByVector(MyArray v);
    void set_diag_elements();
    std::vector<float> get_diag_elements();
    int get_size();
    int get_x(int index);
    int get_y(int index);
    float get_value(int index);
    int* get_x();
    int* get_y();
    float* get_data();
    void ShowData();
    bool isSorted();

    void set_value(int row, int col, float value);
    void pop(int row, int col);
    int CountNonZero();
    SparseMatrixCOO DeleteZeros();
    void DeleteZerosOrder();
	void Show();
    void ShowAsMatrix(int start, int end, int n);
    void ShowAsMatrixNumber(int start, int end, int n);
    void ShowAsMatrixSymbol(int start, int end, int n);
private:
	std::vector<Triplet> v;
    std::vector<float> diag_elems;
	int nonzero;
    int sparse_size;
    bool sorted;

	int *x;
	int *y;
    float *data;
};
int CountNonZero(std::vector<Triplet> t);

#endif
