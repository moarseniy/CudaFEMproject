#include "fem_utils_kernels.h"

// nxm * kxp
__device__
void mult(float *a, float *b, float *tgt,
          size_t n, size_t m,
          size_t k, size_t p,
          bool a_tr, size_t thread_id) {
  float a1;
  if (a_tr) {
    size_t temp = n;
    n = m;
    m = temp;
  }
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < p; j++) {
      tgt[j + i * p + n * p * thread_id] = 0.f;
      for (size_t t = 0; t < n; t++) {
        a1 = a_tr ? a[i + t * m + n * m * thread_id] :
                    a[t + i * m + n * m * thread_id];
        tgt[j + i * p + n * p * thread_id] += a1 * b[j + t * p];
      }
    }
  }
}

__global__
void kernelGenerateMask(float *elements, float *mask, size_t dim, size_t n) {
  size_t index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < n) {

  }
}

void genMask_ker(float *mask, float *elements, size_t dim, size_t elementsCount) {
  kernelGenerateMask<<<(elementsCount + 255) / 256, 256>>>(elements, mask, dim, elementsCount);
}
