
#include <matrix_pack/matrix_pack.h>

Matrix::Matrix(size_t numRows, size_t numCols, float *data) :
  _numRows(numRows),
  _numCols(numCols),
  _numElements(numRows * numCols),
  _data(data) {}

size_t Matrix::get_numRows() const {
  return _numRows;
}

size_t Matrix::get_numCols() const {
  return _numCols;
}

size_t Matrix::get_numElements() const {
  return _numElements;
}

float* Matrix::get_data() const {
  return _data;
}

void Matrix::fillRandValues(float v1, float v2) {
  assert(_numElements > 0);
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      _data[j + i * _numCols] = v1 + rand() / (float) RAND_MAX * (v2 - v1);
    }
  }
}

void Matrix::Show() {
  assert(_numElements > 0);
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      std::cout << _data[j + i * _numCols] << " ";
    }
    std::cout << "\n";
  }
}
