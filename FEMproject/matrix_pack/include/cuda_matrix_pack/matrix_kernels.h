#ifndef MATRIX_KERNELS_H_
#define MATRIX_KERNELS_H_

#include <cstdio>
#include <cuda_runtime.h>
#include <cuda.h>

void multiplyByVec_Ker(const float *data, const size_t numMatr,
                       const float *vec_data, const size_t vecSize,
                       float *tgt);

void copy_Ker(const float *src, float *dst, const size_t numElements);


#endif /* MATRIX_KERNELS_H_ */
