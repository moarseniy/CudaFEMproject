#include <matrix_pack/matrix_pack.h>

CPU_Matrix::CPU_Matrix() :
  Matrix(0, 0, true, nullptr, CPU) {}

CPU_Matrix::CPU_Matrix(size_t numRows, size_t numCols) :
  Matrix(numRows, numCols, true, nullptr, CPU) {
  _data = numRows * numCols > 0 ? new float[_numElements] : nullptr;
}

CPU_Matrix::CPU_Matrix(float *data, size_t numRows, size_t numCols, bool hasData) :
  Matrix(numRows, numCols, hasData, data, CPU) {}

CPU_Matrix::CPU_Matrix(const CPU_Matrix &like) :
  Matrix(like._numRows, like._numCols, true, nullptr, CPU) {
  _data = like._numRows * like._numCols > 0 ? new float[_numElements] : nullptr;
}

CPU_Matrix::~CPU_Matrix() {
  if (_numElements > 0 && _hasData)
    delete[] _data;
}

void CPU_Matrix::product(Matrix &src, Matrix &tgt, bool a_tr, bool b_tr) {
  assert(_numElements > 0 && src.get_numElements() > 0);
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

void CPU_Matrix::product(Matrix &src, bool a_tr, bool b_tr) {
  product(src, *this, a_tr, b_tr);
}

void CPU_Matrix::transpose() {
  std::swap(_numRows, _numCols);
}

float CPU_Matrix::dotProduct(Matrix &vec) {
  assert(_numElements > 0 && _numElements == vec.get_numElements());

  float dotProd = 0.f;
  for (size_t i = 0; i < _numElements; ++i) {
    dotProd += _data[i] * vec.get_data()[i];
  }
  return dotProd;
}

void CPU_Matrix::divideElementwise(float value) {
  assert(_numElements > 0 && _numElements);
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] /= value;
  }
}

void CPU_Matrix::divideElementwise(Matrix &vec, Matrix &tgt) {
  assert(_numElements > 0 && _numElements == vec.get_numElements());
  for (size_t i = 0; i < _numElements; ++i) {
     tgt.get_data()[i] = _data[i] / vec.get_data()[i];
  }
}

void CPU_Matrix::divideElementwise(Matrix &vec) {
  assert(_numElements > 0 && _numElements == vec.get_numElements());
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] /= vec.get_data()[i];
  }
}

void CPU_Matrix::scale(float value) {
  assert(_numElements > 0);
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] *= value;
  }
}

void CPU_Matrix::resize(Matrix &like) {
  this->resize(like.get_numRows(), like.get_numCols());
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

void CPU_Matrix::add(Matrix &src) {
  assert(_numElements > 0 && _numElements == src.get_numElements());
  for (size_t i = 0; i < _numElements; ++i) {
     _data[i] += src.get_data()[i];
  }
}

void CPU_Matrix::sort_by_key(Matrix &keys) {
  assert(_numElements > 0 && _numElements == keys.get_numElements());

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
  assert(_numElements > 0 && _numElements == keys.get_numElements());

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
  assert(_numElements > 0);
  Matrix::copy(target);
  std::stable_sort(target.get_data(),
                   target.get_data() + _numElements);
}

void CPU_Matrix::sort() {
  assert(_numElements > 0);
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
  if (_numRows == _numCols && _numRows == 3) {
    return det3x3(_data);
  } else if (_numRows == _numCols && _numRows == 4) {
    return det4x4(_data);
  } else {
    throw std::runtime_error("CPU_Matrix::det -> No size for det");
  }
}

void CPU_Matrix::setTo(float value) {
  assert(_numElements > 0);
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      _data[j + i * _numCols] = value;
    }
  }
}

void CPU_Matrix::multiplyByVec(const Matrix &vec, Matrix &target) const {
  assert(_numElements > 0);

  size_t numSubMatr = _numElements / _numCols;
  size_t subVec_size = vec.get_numElements() / numSubMatr;
  //  assert(_numCols == vec_size);

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
  assert(_numRows == b.get_numRows() && _numCols == b.get_numCols());
  for (size_t i = 0; i < _numRows; ++i) {
    for (size_t j = 0; j < _numCols; ++j) {
      _data[j + i * _numCols] = alpha * _data[j + i * _numCols] + beta * b(i, j);
    }
  }
}

void CPU_Matrix::copy(Matrix &src, Matrix &tgt) {
  assert(src.get_numElements() > 0);
  if (!src.isSameAs(tgt))
    tgt.resize(src);
  std::memcpy(tgt.get_data(),
              src.get_data(),
              src.get_numElements() * sizeof(float));
}

//float& CPU_Matrix::operator [](size_t index) {
//  return _data[index];
//}
