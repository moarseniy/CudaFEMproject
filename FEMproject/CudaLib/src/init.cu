
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/inner_product.h>
#include <thrust/device_ptr.h>

#include <iostream>
#include <cuda.h>
#include <cusolverSp.h>

#include <math.h>
#include <vector>
#include "init.h"
//#include "femfunc.h"

#include <thrust/reduce.h>

#include "Linal2.h"

using namespace std;

#define CUDA_CHECK_ERROR(err) \
if (err != cudaSuccess) { \
printf("Cuda error: %s\n", cudaGetErrorString(err)); \
printf("Error in file: %s, line: %i\n", __FILE__, __LINE__); \
} \


__global__ void kernelAddWeighted(int n, float *a, float *b, float v1, float v2) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
      a[i] = v1 * a[i] + v2 * b[i];
}

void gpuAddWeighted(float *a, float *b, float v1, float v2, int size) {
  float *d_a, *d_b;

  cudaMalloc(&d_a, size * sizeof(float));
  cudaMalloc(&d_b, size * sizeof(float));

  cudaMemcpy(d_a, a, size * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, size * sizeof(float), cudaMemcpyHostToDevice);

  kernelAddWeighted<<<(size + 255) / 256, 256>>>(size, d_a, d_b, v1, v2);

  cudaMemcpy(a, d_a, size * sizeof(float), cudaMemcpyDeviceToHost);

  cudaFree(d_a);
  cudaFree(d_b);
}

void TEST_THRUST() {
const int N = 7;
int A[N] = {1, 3, 3, 3, 2, 2, 1}; // input keys
int B[N] = {9, 8, 7, 6, 5, 4, 3}; // input values
int C[N];                         // output keys
int D[N];                         // output values
thrust::pair<int*,int*> new_end;
thrust::equal_to<int> binary_pred;
thrust::plus<int> binary_op;
new_end = thrust::reduce_by_key(A, A + N, B, C, D, binary_pred, binary_op);
for (int i = 0; i < 4; ++i) {
    std::cout << C[i] << " ";
}
std::cout << std::endl;
for (int i = 0; i < 4; ++i) {
    std::cout << D[i] << " ";
}
}

void gpuReductionWithMaskAndTransform(float *v, float *mask, int size, float *res, int size_new) {
    thrust::device_vector<float> d_v(v, v + size);
    thrust::device_vector<float> d_mask(mask, mask + size);

    thrust::device_vector<float> d_res(size_new);
    thrust::device_vector<float> d_mask_new(size_new);

    thrust::sort_by_key(d_mask.begin(), d_mask.end(), d_v.begin());
    thrust::reduce_by_key(d_mask.begin(), d_mask.end(), d_v.begin(), d_mask_new.begin(), d_res.begin());

    for (int i = 0; i < size; ++i) {
        float tmp = d_res[static_cast<int>(mask[i])];
        res[i] = tmp;
        //printf("%f ", res[i]);
    }
}

void gpuReductionWithMask(float *v, float *mask, int size, float *res, int size_new) {
    thrust::device_vector<float> d_v(v, v + size);
    thrust::device_vector<float> d_mask(mask, mask + size);

    thrust::device_vector<float> d_res(size_new);
    thrust::device_vector<float> d_mask_new(size_new);

    thrust::sort_by_key(d_mask.begin(), d_mask.end(), d_v.begin());
    thrust::reduce_by_key(d_mask.begin(), d_mask.end(), d_v.begin(), d_mask_new.begin(), d_res.begin());

    for (int i = 0; i < size_new; ++i) {
        float tmp = d_res[i];
        res[i] = tmp;
    }
}

float gpuDotProduct(float *a, float *b, int size) {
    thrust::device_vector<float> d_a(a, a + size);
    thrust::device_vector<float> d_b(b, b + size);
    return thrust::inner_product(thrust::cuda::par, d_a.begin(), d_a.end(), d_b.begin(), 0.0f);
}


