#include <cuda_matrix_pack/cuda_matrix.h>

CUDA_Matrix::CUDA_Matrix() :
  Matrix(0, 0, true, nullptr, CUDA) {}

CUDA_Matrix::CUDA_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, true, nullptr, CUDA) {

  assert(numRows * numCols > 0);
  cudaMalloc(&_data, _numElements * sizeof(float));
}

CUDA_Matrix::CUDA_Matrix(float *data, size_t numRows, size_t numCols, bool hasData) :
  Matrix(numRows, numCols, hasData, data, CUDA) {}

CUDA_Matrix::CUDA_Matrix(const CUDA_Matrix &like) :
  Matrix(like._numRows, like._numCols, true, nullptr, CUDA) {

  assert(like._numRows * like._numCols > 0);
  cudaMalloc(&_data, _numElements * sizeof(float));
}

CUDA_Matrix::~CUDA_Matrix() {
  if (_numElements> 0 && _hasData)
    gpuClearData();
}

void CUDA_Matrix::gpuClearData() {
  cudaFree(_data);
}

void CUDA_Matrix::gpuInitData() {
}

//Thrust-based functions

float CUDA_Matrix::dotProduct(Matrix &vec) {
  assert(_numElements > 0);
  return thrust::inner_product(_data,
                               _data + _numElements,
                               vec.get_data(), 0.f);
}

void CUDA_Matrix::divideElementwise(Matrix &vec) {
  assert(_numElements > 0);
  thrust::transform(_data,
                    _data + _numElements,
                    vec.get_data(), _data,
                    thrust::divides<float>());
}


void CUDA_Matrix::copy(Matrix &src, Matrix &tgt) {
  assert(src.get_numElements() > 0 && src.get_numElements() == tgt.get_numElements());

  if (src.get_device() == CUDA && tgt.get_device() == CUDA) {
    cudaMemcpy(tgt.get_data(), src.get_data(), src.get_numElements() * sizeof(float), cudaMemcpyDeviceToDevice);
  } else if (src.get_device() == CUDA && tgt.get_device() == CPU) {
    cudaMemcpy(tgt.get_data(), src.get_data(), src.get_numElements() * sizeof(float), cudaMemcpyDeviceToHost);
  } else if (src.get_device() == CPU && tgt.get_device() == CUDA) {
    cudaMemcpy(tgt.get_data(), src.get_data(), src.get_numElements() * sizeof(float), cudaMemcpyHostToDevice);
  } else {
    throw std::runtime_error("CUDA_Matrix::copy -> No type for device");
  }
}

void CUDA_Matrix::sort_by_key(Matrix &keys) {
  assert(_numElements > 0);
  thrust::sort_by_key(keys.get_data(),
                      keys.get_data() + keys.get_numElements(),
                      _data);
}

void CUDA_Matrix::reduce_by_key(Matrix &keys, Matrix &target) {
  assert(_numElements > 0);
  thrust::reduce_by_key(keys.get_data(),
                        keys.get_data() + keys.get_numElements(),
                        _data, thrust::make_discard_iterator(),
                        target.get_data());
}

void CUDA_Matrix::sort() {
  assert(_numElements > 0);
  thrust::sort(_data, _data + _numElements);
}

void CUDA_Matrix::setTo(float value) {
  assert(_numElements > 0);
  assert(false);
}

void CUDA_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  assert(_numElements > 0);

  size_t subVec_size = _numCols;
  size_t numSubMatr = vec.get_numElements() / subVec_size;
  multiplyByVec_Ker(_data, numSubMatr, vec.get_data(), subVec_size, target.get_data());
}

void CUDA_Matrix::add(Matrix &src) {
  assert(false);
}

void CUDA_Matrix::addWeighted(Matrix &b, float alpha, float beta) {
  assert(false);
}

void CUDA_Matrix::divideElementwise(float value) {
  assert(false);
}

void CUDA_Matrix::scale(float value) {
  assert(false);
}

float CUDA_Matrix::det() {
  assert(false);
  return 0.f;
}

void CUDA_Matrix::product(Matrix &src, Matrix &tgt) {
  assert(false);
}

void CUDA_Matrix::resize(size_t numRows, size_t numCols) {
  if (numRows != _numRows || numCols != _numCols) {
    if (_numElements != numRows * numCols) {
      if (_numElements > 0) {
        cudaFree(_data);
        _data = NULL;
      }
      if (numRows * numCols > 0) {
        cudaMalloc(&_data, numRows * numCols * sizeof(float));
      } else {
        _data = NULL;
      }
    }

    _numRows = numRows;
    _numCols = numCols;
    _numElements = numRows * numCols;
  }
}
