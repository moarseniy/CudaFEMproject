#include <matrix_pack_cuda/gpu_matrix.h>

CUDA_Matrix::CUDA_Matrix() :
  Matrix(0, 0, nullptr) {}

CUDA_Matrix::CUDA_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, nullptr) {

  assert(numRows * numCols > 0);
  cudaMalloc(&_data, _numElements * sizeof(float));
}

CUDA_Matrix::CUDA_Matrix(const CUDA_Matrix &like) :
  Matrix(like._numRows, like._numCols, nullptr) {

  assert(like._numRows * like._numCols > 0);
  cudaMalloc(&_data, _numElements * sizeof(float));
}

CUDA_Matrix::~CUDA_Matrix() {
  if (_numElements> 0)
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


void CUDA_Matrix::copy(Matrix &target, bool isGPU) const {
  if (isGPU)
    cudaMemcpy(target.get_data(), _data, _numElements * sizeof(float), cudaMemcpyDeviceToHost);
  else
    cudaMemcpy(target.get_data(), _data, _numElements * sizeof(float), cudaMemcpyHostToHost);
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

void CUDA_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  assert(_numElements > 0);

  size_t subVec_size = _numCols;
  size_t numSubMatr = vec.get_numElements() / subVec_size;
  multiplyByVec_Ker(_data, numSubMatr, vec.get_data(), subVec_size, target.get_data());
}

void CUDA_Matrix::product(Matrix &src, Matrix &tgt) {
  assert(false);
}
