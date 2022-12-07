#include <matrix_pack/matrix_pack.h>

CPU_Matrix::CPU_Matrix() :
  Matrix(0, 0, nullptr) {}

CPU_Matrix::CPU_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, nullptr) {
  _data = numRows * numCols > 0 ? new float[_numElements] : nullptr;
}

CPU_Matrix::CPU_Matrix(const CPU_Matrix &like) :
  Matrix(like._numRows, like._numCols, nullptr) {
  _data = like._numRows * like._numCols > 0 ? new float[_numElements] : nullptr;
}

CPU_Matrix::~CPU_Matrix() {
  if (_numElements > 0)
    delete[] _data;
}

void CPU_Matrix::product(Matrix &src, Matrix &tgt) {
  assert(_numElements > 0);
  if (_numCols == src.get_numRows()) {
    for (size_t i = 0; i < _numRows; ++i) {
      for (size_t j = 0; j < src.get_numCols(); ++j) {
        tgt.get_data()[j + tgt.get_numCols() * i] = 0.f;
        for (size_t k = 0; k < _numCols; ++k) {
          tgt.get_data()[j + tgt.get_numCols() * i] += _data[k + i * _numCols] * src.get_data()[j + src.get_numCols() * k];
        }
      }
    }
  } else {
    throw std::runtime_error("product: numCols1 != numRows2\n");
  }
}

//Thrust-based functions

float CPU_Matrix::dotProduct(Matrix &vec) {
  assert(_numElements > 0);
  return thrust::inner_product(thrust::host, _data,
                               _data + _numElements,
                               vec.get_data(), 0.f);
}

void CPU_Matrix::divideElementwise(Matrix &vec) {
  assert(_numElements > 0);
  thrust::transform(thrust::host, _data,
                    _data + _numElements,
                    vec.get_data(), _data,
                    thrust::divides<float>());
}

void CPU_Matrix::sort_by_key(Matrix &keys) {
  assert(_numElements > 0);
  thrust::sort_by_key(thrust::host, keys.get_data(),
                      keys.get_data() + keys.get_numElements(), _data);
}

void CPU_Matrix::reduce_by_key(Matrix &keys, Matrix &target) {
  assert(_numElements > 0);
  thrust::reduce_by_key(thrust::host, keys.get_data(),
                        keys.get_data() + keys.get_numElements(),
                        _data, thrust::make_discard_iterator(),
                        target.get_data());
}

void CPU_Matrix::sort() {
  assert(_numElements > 0);
  std::stable_sort(_data, _data + _numElements);
}

void CPU_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  assert(_numElements > 0);

  size_t subVec_size = _numCols;
  size_t numSubMatr = vec.get_numElements() / subVec_size;
//  assert(_numCols == vec_size);

  for (size_t index = 0; index < numSubMatr; ++index) {
    for (size_t j = 0; j < subVec_size; ++j) {
      target.get_data()[j + index * subVec_size] = 0.f;
      for (size_t k = 0; k < subVec_size; ++k) {
        target.get_data()[j + index * subVec_size] += _data[k + j * subVec_size + index * subVec_size * subVec_size] *
           vec.get_data()[k + index * subVec_size];
      }
    }
  }
}

void CPU_Matrix::copy(Matrix &target, bool isGPU) const {
  if (isGPU)
    cudaMemcpy(target.get_data(), _data, _numElements * sizeof(float), cudaMemcpyHostToDevice);
  else
    cudaMemcpy(target.get_data(), _data, _numElements * sizeof(float), cudaMemcpyHostToHost);
}
