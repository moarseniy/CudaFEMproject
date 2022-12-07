#include "matrix_kernels.h"

__global__ void kernelMultiply(const size_t n, const float *data, const float *vec, const size_t size, float *tgt) {
  size_t index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < n) {
    for (size_t j = 0; j < size; ++j) {
      tgt[j + index * size] = 0.f;
      for (size_t k = 0; k < size; ++k) {
        tgt[j + index * size] += data[k + j * size + index * size * size] * vec[k + index * size];
      }
    }
  }
}

void multiplyByVec_Ker(const float *data, const size_t numMatr,
                       const float *vec_data, const size_t vecSize,
                       float *tgt) {
  kernelMultiply<<<(numMatr + 255) / 256, 256>>>(numMatr, data, vec_data, vecSize, tgt);
}
