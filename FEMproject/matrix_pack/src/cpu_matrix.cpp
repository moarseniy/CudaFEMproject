#include <matrix_pack/matrix_pack.h>

CPU_Matrix::CPU_Matrix() :
  Matrix(0, 0, true, nullptr, CPU) {}

CPU_Matrix::CPU_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, true, nullptr, CPU) {
  _data = numRows * numCols > 0 ? new float[_numElements] : nullptr;
}

CPU_Matrix::CPU_Matrix(float *data, size_t numRows, size_t numCols, bool hasData) :
  Matrix(numRows, numCols, hasData, data, CPU) {}

CPU_Matrix::CPU_Matrix(const CPU_Matrix &like, bool copy) :
  Matrix(like._numRows, like._numCols, true, nullptr, CPU) {
  _data = like._numRows * like._numCols > 0 ? new float[_numElements] : nullptr;
  if (copy) {
    CPU_Matrix::copy(like, *this);
  }
}

CPU_Matrix::~CPU_Matrix() {
  if (_numElements > 0 && _hasData)
    delete[] _data;
}

void CPU_Matrix::product(Matrix &src, Matrix &tgt, bool a_tr, float scaleA,
                                                   bool b_tr, float scaleB) {
  CheckAssert(_numElements > 0 && src.get_numElements() > 0);
  if (a_tr)
    this->transpose();
  if (b_tr)
    src.transpose();

  float a, b;
  if (_numCols == src.get_numRows()) {
    for (size_t i = 0; i < _numRows; ++i) {
      for (size_t j = 0; j < src.get_numCols(); ++j) {
        tgt.get_data()[j + tgt.get_numCols() * i] = 0.f;
        for (size_t k = 0; k < _numCols; ++k) {
          a = a_tr ? _data[i + k * _numRows] :
                     _data[k + i * _numCols];
          b = b_tr ? src.get_data()[j + src.get_numCols() * k] :
                     src.get_data()[j + src.get_numCols() * k]; // TODO: FIX IT!
          tgt.get_data()[j + tgt.get_numCols() * i] += a * b;
        }
      }
    }
  } else {
    throw std::runtime_error("CPU_Matrix::product: numCols1 != numRows2\n");
  }

  if (a_tr)
    this->transpose();
  if (b_tr)
    src.transpose();
}

void CPU_Matrix::bmm(const Matrix &a, size_t subRowsA, size_t subColsA, bool trans_a,
                     const Matrix &b, size_t subColsB, bool trans_b, size_t batchCount, const float alpha) {
  CheckAssert(a.get_numElements() > 0 && b.get_numElements() > 0 && batchCount > 0);
  resize(batchCount, subRowsA * subColsB);
  for (size_t p = 0; p < batchCount; ++p) {
    size_t strideA = subRowsA * subColsA == a.get_numElements() ? 0 : p * subRowsA * subColsA;
    size_t strideB = subColsA * subColsB == b.get_numElements() ? 0 : p * subColsA * subColsB;
    size_t strideC = p * subRowsA * subColsB;
    for (size_t m = 0; m < subRowsA; ++m) {
      for (size_t n = 0; n < subColsB; ++n) {
        size_t c_id = m + n * subRowsA + strideC;
        float c_temp = 0.f;
        for (size_t k = 0; k < subColsA; ++k) {
          size_t a_id = trans_a ? k + m * subColsA + strideA : m + k * subRowsA + strideA;
          size_t b_id = trans_b ? n + k * subColsB + strideB : k + n * subColsA + strideB;
          c_temp += a.get_data()[a_id] * b.get_data()[b_id];
        }
        _data[c_id] = alpha * c_temp;
      }
    }
  }
}

void CPU_Matrix::transpose() {
  std::swap(_numRows, _numCols);
}

float CPU_Matrix::dotProduct(Matrix &vec) {
  CheckAssert(_numElements > 0 && _numElements == vec.get_numElements());

  float dotProd = 0.f;
  for (size_t i = 0; i < _numElements; ++i) {
    dotProd += _data[i] * vec.get_data()[i];
  }
  return dotProd;
}

void CPU_Matrix::divideElementwise(float value) {
  CheckAssert(_numElements > 0);
  this->scale(1.f / value); // TODO: Change to thrust or cublas
}

void CPU_Matrix::divideElementwise(Matrix &vec, Matrix &tgt) {
  CheckAssert(_numElements > 0 && _numElements == vec.get_numElements());
  for (size_t i = 0; i < _numElements; ++i) {
     tgt.get_data()[i] = _data[i] / vec.get_data()[i];
  }
}

void CPU_Matrix::divideElementwise(Matrix &vec) {
  CheckAssert(_numElements > 0 && _numElements == vec.get_numElements());
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] /= vec.get_data()[i];
  }
}

void CPU_Matrix::divideElementwise(Matrix &v, Axis ax) {
  CheckAssert(_numElements > 0 && v.get_numElements() > 0);
  CPU_Matrix inverted(dynamic_cast<CPU_Matrix&>(v));
  inverted.setTo(1.f);
  inverted.divideElementwise(v);
  scale(inverted, ax);
}

void CPU_Matrix::scale(float value) {
  CheckAssert(_numElements > 0);
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] *= value;
  }
}

void CPU_Matrix::scale(Matrix &vec, Axis ax) {
  CheckAssert(_numElements > 0 && vec.get_numElements() > 0);
  switch (ax) {
  case X:
    CheckAssert(_numRows == vec.get_numElements());
    for (size_t i = 0; i < _numRows; ++i) {
      for (size_t j = 0; j < _numCols; ++j) {
          _data[j + i * _numCols] *= vec[i];
      }
    }
    break;
  case Y:
    CheckAssert(_numCols == vec.get_numElements());
    for (size_t i = 0; i < _numRows; ++i) {
      for (size_t j = 0; j < _numCols; ++j) {
          _data[j + i * _numCols] *= vec[j];
      }
    }
    break;
  case ALL:
    CheckAssert(_numElements == vec.get_numElements());
    for (size_t i = 0; i < _numElements; ++i) {
      _data[i] *= vec[i];
    }
    break;
  }
}

void CPU_Matrix::resize(const Matrix &like) {
  resize(like.get_numRows(), like.get_numCols());
}

void CPU_Matrix::resize(size_t numRows, size_t numCols) {
  if (numRows != _numRows || numCols != _numCols) {
    if (_numElements != numRows * numCols) {
      delete[] _data;
      _data = new float[numRows * numCols];
    }
    _numRows = numRows;
    _numCols = numCols;
    _numElements = numRows * numCols;
  }
}

void CPU_Matrix::flatten() {
  this->resize(1, _numElements);
}

void CPU_Matrix::add(Matrix &src) {
  CheckAssert(_numElements > 0 && _numElements == src.get_numElements());
  this->addWeighted(src, 1.f, 1.f);
}

void CPU_Matrix::sort_by_key(Matrix &keys) {
  CheckAssert(_numElements > 0 && _numElements == keys.get_numElements());

  std::multimap<float, float> multimap_data;
  for (int i = 0; i < _numElements; ++i) {
    multimap_data.insert(std::multimap<float, float>::value_type(keys.get_data()[i], _data[i]));
  }

  size_t i = 0;
  for (auto &el : multimap_data) {
    keys.get_data()[i] = el.first;
    _data[i] = el.second;
    ++i;
  }
}

void CPU_Matrix::reduce_by_key(Matrix &keys, Matrix &target) {
  CheckAssert(_numElements > 0 && _numElements == keys.get_numElements());

  std::multimap<float, float> multimap_data;
  std::set<float> unique_keys;
  for (int i = 0; i < _numElements; ++i) {
    multimap_data.insert(std::multimap<float, float>::value_type(keys.get_data()[i], _data[i]));
    unique_keys.insert(keys.get_data()[i]);
  }

  size_t i = 0;
  for (auto &key : unique_keys) {
    auto range = multimap_data.equal_range(key);
    target.get_data()[i] = std::accumulate(range.first, range.second,
        0.f,
        [](float a, std::pair<float, float> b) { return a + b.second; });
    ++i;
  }
}

void CPU_Matrix::sort(Matrix &target) {
  CheckAssert(_numElements > 0);
  Matrix::copy(target);
  std::stable_sort(target.get_data(),
                   target.get_data() + target.get_numElements());
}

void CPU_Matrix::sort() {
  CheckAssert(_numElements > 0);
  std::stable_sort(_data, _data + _numElements);
}

float det3(float a0, float a1, float a2,
             float a3, float a4, float a5,
             float a6, float a7, float a8) {
  return a0 * a4 * a8 +
      a1 * a6 * a5 +
      a2 * a3 * a7 -
      a6 * a4 * a2 -
      a0 * a5 * a7 -
      a1 * a3 * a8;
}

float det3x3(float *c) {
  return c[0] * c[4] * c[8] +
        c[1] * c[6] * c[5] +
        c[2] * c[3] * c[7] -
        c[6] * c[4] * c[2] -
        c[0] * c[5] * c[7] -
        c[1] * c[3] * c[8];
}

float det4x4(float *m) {
  float v1 = det3(m[5], m[6], m[7], m[9],
      m[10], m[11], m[13], m[14], m[15]);
  float v2 = det3(m[1], m[2], m[3], m[9],
      m[10], m[11], m[13], m[14], m[15]);
  float v3 = det3(m[1], m[2], m[3], m[5],
      m[6], m[7], m[13], m[14], m[15]);
  float v4 = det3(m[1], m[2], m[3], m[5],
      m[6], m[7], m[9], m[10], m[11]);
  return v1 - v2 + v3 - v4;
}

float CPU_Matrix::det() {
  CheckAssert(_numRows == _numCols);
  if (_numRows == 3) {
    return det3x3(_data);
  } else if (_numRows == 4) {
    return det4x4(_data);
  } else {
    throw std::runtime_error("CPU_Matrix::det -> No size for det");
  }
}

void CPU_Matrix::setTo(float value) {
  CheckAssert(_numElements > 0);
  for (size_t i = 0; i < _numElements; ++i) {
    _data[i] = value;
  }
}

void CPU_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  CheckAssert(_numElements > 0);

  size_t numSubMatr = _numElements / _numCols;
  size_t subVec_size = vec.get_numElements() / numSubMatr;
  //  CheckAssert(_numCols == vec_size);

  for (size_t index = 0; index < numSubMatr; ++index) {
    for (size_t j = 0; j < subVec_size; ++j) {
      target.get_data()[j + index * subVec_size] = 0.f;
      for (size_t k = 0; k < subVec_size; ++k) {
        target.get_data()[j + index * subVec_size] +=
           _data[k + j * subVec_size + index * subVec_size * subVec_size] *
           vec.get_data()[k + index * subVec_size];
      }
    }
  }
}

void CPU_Matrix::addWeighted(Matrix &b, float alpha, float beta) {
  CheckAssert(_numRows == b.get_numRows() && _numCols == b.get_numCols());
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      _data[j + i * _numCols] = alpha * _data[j + i * _numCols] + beta * b(i, j);
    }
  }
}

void CPU_Matrix::copy(Matrix &tgt) {
  CPU_Matrix::copy(dynamic_cast<Matrix&>(*this), tgt);
}

void CPU_Matrix::copy(const Matrix &src, Matrix &tgt) {
  CheckAssert(src.get_numElements() > 0);
  if (!src.isSameAs(tgt))
    tgt.resize(src);
  if (src.get_device() == CPU && tgt.get_device() == CPU) {
    std::memcpy(tgt.get_data(),
                src.get_data(),
                src.get_numElements() * sizeof(float));
  } else {
    src.copy(tgt);
  }
}

//float& CPU_Matrix::operator [](size_t index) {
//  return _data[index];
//}
