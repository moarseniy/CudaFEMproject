#ifndef MATRIX_KERNELS_H_
#define MATRIX_KERNELS_H_

#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <curand.h>
#include <curand_kernel.h>

#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/inner_product.h>
#include <thrust/fill.h>

float thrust_dotProduct_Ker(float *src1, float *src2, size_t size);
void thrust_divideElementwise_Ker(float *src1, float *src2, float *tgt, size_t size);
void thrust_sortByKey_Ker(float *keys, float *src, size_t size);
void thrust_reduceByKey_Ker(float *keys, float *src, float *tgt, size_t size);
void thrust_sort_Ker(float *src, size_t size);
void thrust_setTo_Ker(float *src, size_t size, float v);

void multiplyByVec_Ker(const size_t numMatr, const float *data,
                       const float *vec_data, const size_t vecSize,
                       float *tgt);

void setTo_Ker(const size_t size, float *data, const float value);
void addWeighted_Ker(const size_t size, float *data, const float *src,
                     const float alpha, const float beta);
void scale_Ker(const size_t size, float *data, const float value);
void runTest();
void runTest(float *a, float *b, float *c);
void runTest2();
void runTest3();
void runTest4();
#endif /* MATRIX_KERNELS_H_ */
