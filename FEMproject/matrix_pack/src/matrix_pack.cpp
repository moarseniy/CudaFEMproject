
#include <matrix_pack/matrix_pack.h>
#ifdef WITH_CUDA
#include <cuda_matrix_pack/cuda_matrix.h>
#endif

Matrix::Matrix(size_t numRows, size_t numCols, bool hasData, float *data, DEVICE_NAME device) :
  _device(device),
  _numRows(numRows),
  _numCols(numCols),
  _numElements(numRows * numCols),
  _data(data),
  _hasData(hasData) {}

DEVICE_NAME Matrix::get_device() const {
  return _device;
}

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

void Matrix::copy(Matrix &tgt) const {
  if (this->get_device() == CPU &&
      tgt.get_device() == CPU) {
    CPU_Matrix::copy(*this, tgt);
#ifdef WITH_CUDA
  } else if (this->get_device() == CUDA ||
             tgt.get_device() == CUDA) {
    CUDA_Matrix::copy(*this, tgt);
#elif WITH_ONEAPI
  } else if (this->get_device() == ONEAPI ||
             tgt.get_device() == ONEAPI) {
    assert(false);
#endif
  } else {
    throw std::runtime_error("Matrix::copy -> No type for device");
  }
}

void Matrix::setDevice(DEVICE_NAME device) {
  _device = device;
}

void Matrix::Show() {
  assert(_numElements > 0);
//  assert(_device == CPU);
  CPU_Matrix tmp;
  this->copy(tmp);
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      std::cout << tmp.get_data()[j + i * _numCols] << " ";
    }
    std::cout << "\n";
  }
}

bool Matrix::isSameAs(Matrix &src) const {
  return _numRows == src.get_numRows() && _numCols == src.get_numCols();
}

Matrix* Matrix::setMatrix(DEVICE_NAME device, size_t numRows, size_t numCols) {
  if (device == CPU) {
    return new CPU_Matrix(numRows, numCols);
#ifdef WITH_CUDA
  } else if (device == CUDA) {
    return new CUDA_Matrix(numRows, numCols);
#endif
  } else {
    throw std::runtime_error("Matrix::setMatrix -> No type for device");
  }
}

Matrix* Matrix::setMatrix(DEVICE_NAME device) {
  return setMatrix(device, 0, 0);
}

Matrix* Matrix::setMatrix(Matrix &like) {
  return setMatrix(like.get_device(), like.get_numRows(), like.get_numCols());
}

Matrix* Matrix::setVector(DEVICE_NAME device, size_t size) {
  return setMatrix(device, size, 1);
}

Matrix* Matrix::setVector(DEVICE_NAME device) {
  return setMatrix(device, 0, 0);
}

Matrix* Matrix::setVector(Matrix &like) {
  // TODO: think about it
  return setVector(like.get_device(), like.get_numCols());
}

std::unique_ptr<Matrix> Matrix::subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol) const {
  if (_device == CPU) {
    return std::unique_ptr<CPU_Matrix>{new CPU_Matrix(_data + startRow * _numCols + startCol, endRow - startRow, endCol - startCol)};
#ifdef WITH_CUDA
  } else if (_device == CUDA) {
    return std::unique_ptr<CUDA_Matrix>{new CUDA_Matrix(_data + startRow * _numCols + startCol, endRow - startRow, endCol - startCol)};
#endif
  } else {
    throw std::runtime_error("Matrix::subMatrix -> No type for device");
  }
}

void Matrix::subMatrix(size_t startRow, size_t endRow, size_t startCol, size_t endCol, Matrix &target) const {
  target.resize(endRow - startRow, endCol - startCol);
  // TODO: complete this case
//  this->copy(target, startRow, endRow, startCol, endCol, 0, 0);
}

std::unique_ptr<Matrix> Matrix::getRows(size_t startRow, size_t endRow) const {
  return subMatrix(startRow, endRow, 0, _numCols);
}

void Matrix::getRow(size_t numRow, Matrix &target) const {
  subMatrix(numRow, numRow + 1, 0, _numCols, target);
}

std::unique_ptr<Matrix> Matrix::getRow(size_t numRow) const {
  return subMatrix(numRow, numRow + 1, 0, _numCols);
}

void Matrix::getRows(size_t startRow, size_t endRow, Matrix &target) const {
  subMatrix(startRow, endRow, 0, _numCols, target);
}

float& Matrix::operator [](size_t index) {
  assert(_device == CPU);
  return _data[index];
}

float& Matrix::operator ()(size_t i, size_t j) {
  assert(_device == CPU);
  return _data[j + i * _numCols];
}
