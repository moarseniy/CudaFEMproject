#include <cuda_matrix.h>

cublasHandle_t CUDA_Matrix::_handle;

void CUDA_Matrix::_initCUDA() {
  if (!_handle)
    cublasCreate(&_handle);
}

CUDA_Matrix::CUDA_Matrix() :
  Matrix(0, 0, true, nullptr, CUDA) {
  _gpuInitData();
}

CUDA_Matrix::CUDA_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, true, nullptr, CUDA) {
  _gpuInitData();
}

CUDA_Matrix::CUDA_Matrix(float *data, size_t numRows, size_t numCols, bool hasData) :
  Matrix(numRows, numCols, hasData, data, CUDA) {
  _gpuInitData();
}

CUDA_Matrix::CUDA_Matrix(const CUDA_Matrix &like, bool copy) :
  Matrix(like._numRows, like._numCols, true, nullptr, CUDA) {
  _gpuInitData();
  if (copy) {
//    like.copy(*this);
    CUDA_Matrix::copy(like, *this);
  }
}

CUDA_Matrix::~CUDA_Matrix() {
  if (_numElements> 0 && _hasData)
    _gpuClearData();
}

void CUDA_Matrix::_gpuClearData() {
  cudaFree(_data);
}

void CUDA_Matrix::_gpuInitData() {
  if (_numElements > 0) {
    cudaMalloc((void **)&_data, _numElements * sizeof(float));
  }
}

//Thrust-based functions

float CUDA_Matrix::dotProduct(Matrix &vec) {
  CheckAssert(_numElements > 0);
  return thrust_dotProduct_Ker(_data, vec.get_data(), _numElements);
}

void CUDA_Matrix::divideElementwise(Matrix &v, Axis ax) {
  CheckAssert(_numElements > 0 && v.get_numElements() > 0);
  CUDA_Matrix inverted(dynamic_cast<CUDA_Matrix&>(v));
  inverted.setTo(1.f);
  inverted.divideElementwise(v);
  scale(inverted, ax);
}

void CUDA_Matrix::divideElementwise(Matrix &vec, Matrix &tgt) {
  CheckAssert(_numElements > 0);
  CheckAssert(_numElements == vec.get_numElements());
  CheckAssert(_numElements == tgt.get_numElements());
  thrust_divideElementwise_Ker(_data, vec.get_data(), tgt.get_data(), _numElements);
}

void CUDA_Matrix::divideElementwise(Matrix &vec) {
  divideElementwise(vec, *this);
}

void CUDA_Matrix::copy(Matrix &tgt) {
  CUDA_Matrix::copy(dynamic_cast<Matrix&>(*this), tgt);
}

void CUDA_Matrix::copy(const Matrix &src, Matrix &tgt) {
  CheckAssert(src.get_numElements() > 0);
  if (!src.isSameAs(tgt))
    tgt.resize(src);

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
  CheckAssert(_numElements > 0);
  CheckAssert(_numElements == keys.get_numElements());
  thrust_sortByKey_Ker(keys.get_data(), _data, _numElements);
}

void CUDA_Matrix::reduce_by_key(Matrix &keys, Matrix &target) {
  CheckAssert(_numElements > 0);
  thrust_reduceByKey_Ker(keys.get_data(), _data, target.get_data(), _numElements);
}

void CUDA_Matrix::sort() {
  CheckAssert(_numElements > 0);
  thrust_sort_Ker(_data, _numElements);
}

void CUDA_Matrix::sort(Matrix &target) {
  CheckAssert(_numElements > 0);
  copy(target);
  thrust_sort_Ker(target.get_data(), target.get_numElements());
}

float CUDA_Matrix::min() {
  CheckAssert(_numElements > 0);
  return thrust_min_Ker(_data, _numElements);
}

float CUDA_Matrix::max() {
  CheckAssert(_numElements > 0);
  return thrust_max_Ker(_data, _numElements);
}

void CUDA_Matrix::setTo(float value) {
  CheckAssert(_numElements > 0); // TODO: compare kernel and thrust
//  setTo_Ker(_numElements, _data, value);
  thrust_setTo_Ker(_data, _numElements, value);
}

void CUDA_Matrix::getDiagonal(size_t size, Matrix &tgt) {
  CheckAssert(_numElements > 0);
  tgt.resize(_numRows, size);

  getDiagonal_Ker(_numRows, size, _data, tgt.get_data());
}

void CUDA_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  CheckAssert(_numElements > 0);

  size_t subVec_size = _numCols;
  size_t numSubMatr = vec.get_numElements() / subVec_size;
  multiplyByVec_Ker(numSubMatr, _data, vec.get_data(), subVec_size, target.get_data());
}

void CUDA_Matrix::add(Matrix &src) {
  addWeighted(src, 1.f, 1.f);
}

void CUDA_Matrix::add(Matrix &src, Matrix &tgt) {
  addWeighted(src, 1.f, 1.f, tgt);
}

void CUDA_Matrix::subtract(Matrix &src) {
  addWeighted(src, 1.f, -1.f);
}

void CUDA_Matrix::subtract(Matrix &src, Matrix &tgt) {
  addWeighted(src, 1.f, -1.f, tgt);
}

void CUDA_Matrix::addWeighted(Matrix &b, float alpha, float beta) {
  addWeighted(b, alpha, beta, *this);
}

void CUDA_Matrix::addWeighted(Matrix &b, float alpha, float beta, Matrix &target) {
  CheckAssert(_numElements > 0 && _numElements == b.get_numElements());
  CheckAssert(this->isSameAs(b));
  if (!this->isSameAs(target))
    target.resize(*this);

  cublasSgeam(_handle, CUBLAS_OP_N, CUBLAS_OP_N,
              _numRows, _numCols,
              &alpha, _data, _numRows,
              &beta, b.get_data(), _numRows,
              target.get_data(), _numRows);
}

void CUDA_Matrix::divideElementwise(float value) {
  CheckAssert(_numElements > 0);
  this->scale(1.f / value);
}

void CUDA_Matrix::scale(float value) {
  CheckAssert(_numElements > 0);
  scale_Ker(_numElements, _data, value);
}

void CUDA_Matrix::scale(Matrix &v, Axis ax) {
  CheckAssert(_numElements > 0 && v.get_numElements() > 0);
  switch (ax) {
  case X:
    CheckAssert(_numRows == v.get_numElements());
    cublasSdgmm(_handle, CUBLAS_SIDE_RIGHT,
                _numCols, _numRows,
                _data, _numCols,
                v.get_data(), 1,
                _data, _numCols);
    break;
  case Y:
    CheckAssert(_numCols == v.get_numElements());
    cublasSdgmm(_handle, CUBLAS_SIDE_LEFT,
                _numCols, _numRows,
                _data, _numCols,
                v.get_data(), 1,
                _data, _numCols);
    break;
  default:
    CheckAssert(_numElements == v.get_numElements());
    cublasSdgmm(_handle, CUBLAS_SIDE_RIGHT,
                1, _numElements,
                _data, 1,
                v.get_data(), 1,
                _data, 1);
    break;
  }
}

float CUDA_Matrix::det() {
  CheckAssert(false);
  return 0.f;
}

float CUDA_Matrix::l2norm() {
  CheckAssert(_numElements > 0);
  float l2norm_result = 0.0f;
  cublasSnrm2(_handle,
              _numElements, _data, 1,
              &l2norm_result);
  return l2norm_result;
}

//void CUDA_Matrix::product(Matrix &src, Matrix &tgt, bool a_tr, bool b_tr) {
//  CheckAssert(false);
//}

void CUDA_Matrix::product(Matrix &src, Matrix &tgt, bool a_trans, float scaleA,
                                                    bool b_trans, float scaleB)  {
//  if (scaleThis == 0.f)
//    this->resize(b.getNumRows(), a.getNumCols());
//  else
//    my_CheckAssert(this->getNumCols() == a.getNumCols() && this->getNumRows() == b.getNumRows());


//  cublasSgemm('n', 'n', _numCols, b.get_numRows(), a.get_numRows(), scaleA,
//              a.get_data(), a.get_numCols(), b.get_data(), b.get_numCols(), scaleB,
//              _data, _numCols);
}

void CUDA_Matrix::bmm(const Matrix &a, size_t subRowsA, size_t subColsA, bool trans_a,
                      const Matrix &b, size_t subColsB, bool trans_b, size_t batchCount, const float alpha) {
  CheckAssert(a.get_numElements() > 0 && b.get_numElements() > 0 && batchCount > 0);
//  resize(batchCount, subRowsA * subColsB);

//  cublasSetMathMode(_handle, CUBLAS_PEDANTIC_MATH);
//  cublasSetComp
  const float beta = 0.f;
  int strideA = subRowsA * subColsA == a.get_numElements() ? 0 : subRowsA * subColsA;
  int strideB = subColsA * subColsB == b.get_numElements() ? 0 : subColsA * subColsB;
  int strideC = subRowsA * subColsB;
  cublasSgemmStridedBatched(_handle,
                            trans_a ? CUBLAS_OP_T : CUBLAS_OP_N,
                            trans_b ? CUBLAS_OP_T : CUBLAS_OP_N,
                            subRowsA, subColsB, subColsA,
                            &alpha,
                            a.get_data(), trans_a ? subColsA : subRowsA, strideA,
                            b.get_data(), trans_b ? subColsB : subColsA, strideB,
                            &beta,
                            _data, subRowsA, strideC,
                            batchCount);

  cudaDeviceSynchronize();
}

void CUDA_Matrix::resize(const Matrix &like) {
  resize(like.get_numRows(), like.get_numCols());
}

void CUDA_Matrix::resize(size_t numRows, size_t numCols) {
  if (numRows != _numRows || numCols != _numCols) {
    if (_numElements != numRows * numCols) {
      if (_numElements > 0) {
        cudaFree(_data);
        _data = NULL;
      }
      if (numRows * numCols > 0) {
        cudaMalloc((void **)&_data, numRows * numCols * sizeof(float));
      } else {
        _data = NULL;
      }
    }

    _numRows = numRows;
    _numCols = numCols;
    _numElements = numRows * numCols;
  }
}

void CUDA_Matrix::flatten() {
  resize(1, _numElements);
}

void CUDA_Matrix::uniformRandomize(float v1, float v2) {
  CheckAssert(_numElements > 0);
  uniformRandomize_Ker(_data, _numElements, v1, v2);
}

void CUDA_Matrix::fillSequence(float startValue) {
  CheckAssert(_numElements > 0);
  fillSequence_ker(_data, _numElements, startValue);
}
