#include <init.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/inner_product.h>
#include <thrust/device_ptr.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/reduce.h>

#include <iostream>
#include <cuda.h>
#include <cusolverSp.h>

#include <math.h>
#include <vector>

#include <matrix_pack/matrix_pack.h>

using namespace std;

#define CUDA_CHECK_ERROR(err) \
  if (err != cudaSuccess) { \
  printf("Cuda error: %s\n", cudaGetErrorString(err)); \
  printf("Error in file: %s, line: %i\n", __FILE__, __LINE__); \
} \


gpuDataKeeper::gpuDataKeeper(int DIM, int elementsCount, int nodesCount, bool doAssemblyRes) : diag(6 * (DIM - 1) * elementsCount), r(6 * (DIM - 1) * elementsCount),
                                                      m(6 * (DIM - 1) * elementsCount), z(6 * (DIM - 1) * elementsCount), s(6 * (DIM - 1) * elementsCount),
                                                      p(6 * (DIM - 1) * elementsCount), u(6 * (DIM - 1) * elementsCount),
                                                      x(6 * (DIM - 1) * elementsCount, 0.0f), mask(6 * (DIM - 1) * elementsCount),
                                                      n_adjelem(DIM * nodesCount), gpuB(3 * (DIM - 1) * 6 * (DIM - 1) * elementsCount, 0.0f), gpuElements((DIM + 1) * elementsCount),
                                                      gpuKlocals(6 * (DIM - 1) * 6 * (DIM - 1) * elementsCount, 0.0f), gpuFlocals(6 * (DIM - 1) * elementsCount, 0.0f), tmp(6 * (DIM - 1) * elementsCount),
                                                      loads(6 * (DIM - 1) * elementsCount, 0.0f)
{
  CheckRunTime(__func__)
  if (doAssemblyRes)
    temp_res.resize(DIM * nodesCount);
}

gpuDataKeeper::~gpuDataKeeper() {}

//WeightedAddCoef::WeightedAddCoef(float v1, float v2): val1(v1), val2(v2) {}
//WeightedAddCoef::WeightedAddCoef(float v1, float v2, float v3): val1(v1), val2(v2), val3(v3) {}

void gpuDataKeeper::setZeroVec() {
  CheckRunTime(__func__)
  thrust::fill(x.begin(), x.end(), 0.0f);
}

gpuDataKeeper_DYN::gpuDataKeeper_DYN(int DIM, int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped) :
  gpuDataKeeper(DIM, elementsCount, nodesCount, doAssemblyRes ),
  vel(3 * DIM * elementsCount, 0.0f), displ(3 * DIM * elementsCount, 0.0f), displ_global(DIM * nodesCount) {
  // default: explicit scheme
  this->SLAU_matrix_coefs.val1 = 1.0f;
  this->SLAU_matrix_coefs.val2 = 0.0f;
  this->SLAU_matrix_coefs.val3 = 0.0f;
  this->isLumped = isLumped;
  if (isLumped) {
    diagM.resize(3 * DIM * elementsCount);
//    gpuMlocals.resize(6 * 6 * elementsCount);
//    thrust::fill(gpuMlocals.begin(), gpuMlocals.end(), 0.0f); // Initialize
  } else {
    gpuMlocals.resize(6 * 6 * elementsCount);
    thrust::fill(gpuMlocals.begin(), gpuMlocals.end(), 0.0f); // Initialize
  }

  // ToDO: initialize displ, vel, x to 0
}

void gpuDataKeeper::copyBmatrixToHost(float *all_B) {
  gpuCopyDeviceToHost(this->get_B(), all_B, this->gpuB.size());
}

gpuDataKeeper_DYN::~gpuDataKeeper_DYN() {}

void gpuDataKeeper_DYN::set_SLAU_matrix_coefs(float cM, float cK) {
  this->SLAU_matrix_coefs.val1 = cM;
  this->SLAU_matrix_coefs.val2 = cK;
}

void gpuDataKeeper_DYN_DAMP::set_SLAU_matrix_coefs(float cM, float cK, float cC) {
  this->SLAU_matrix_coefs.val1 = cM;
  this->SLAU_matrix_coefs.val2 = cK;
  this->SLAU_matrix_coefs.val3 = cC;
}

void gpuDataKeeper_DYN_DAMP::set_damping_coefs(float cM, float cK) {
  this->damping_coefs.val1 = cM;
  this->damping_coefs.val2 = cK;
}

gpuDataKeeper_DYN_DAMP::gpuDataKeeper_DYN_DAMP(int DIM, int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped, float damping_alpha, float damping_beta) :
  gpuDataKeeper_DYN(DIM, elementsCount, nodesCount, doAssemblyRes, isLumped ) {
  this->damping_coefs.val1 = damping_alpha;
  this->damping_coefs.val2 = damping_beta;
}

gpuDataKeeper_DYN_DAMP::~gpuDataKeeper_DYN_DAMP() {}

__global__ void kernelAddWeighted(int n, float *a, float *b, float v1, float v2) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    a[i] = v1 * a[i] + v2 * b[i];
}

__global__ void kernelAdd(int n, float *a, float *b) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    a[i] = a[i] + b[i];
}

void gpuAddWeighted(float *a, float *b, float v1, float v2, int size) {
  CheckRunTime(__func__)
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

void gpuAddWeighted2(float *a_ptr, float *b_ptr,
                     float v1, float v2, int size) {
  CheckRunTime(__func__)
  //  float* a_ptr = thrust::raw_pointer_cast( a.data() );
  //  float* b_ptr = thrust::raw_pointer_cast( b.data() );
  kernelAddWeighted<<<(size + 255) / 256, 256>>>(size, a_ptr, b_ptr, v1, v2);
}

void gpuAdd(float *a_ptr, float *b_ptr, int size) {
  CheckRunTime(__func__)
  //  float* a_ptr = thrust::raw_pointer_cast( a.data() );
  //  float* b_ptr = thrust::raw_pointer_cast( b.data() );
  kernelAdd<<<(size + 255) / 256, 256>>>(size, a_ptr, b_ptr);
}

void gpuDivide(float *a_ptr, float *b_ptr, int size) {
  CheckRunTime(__func__)

  thrust::transform(thrust::device,
                    thrust::device_pointer_cast(a_ptr),
                    thrust::device_pointer_cast(a_ptr + size),
                    thrust::device_pointer_cast(b_ptr),
                    thrust::device_pointer_cast(a_ptr),
                    thrust::divides<float>());
}

void gpuDivide_res(float *a_ptr, float *b_ptr, int size, float *res) {
  CheckRunTime(__func__)

  thrust::transform(thrust::device,
                    thrust::device_pointer_cast(a_ptr),
                    thrust::device_pointer_cast(a_ptr + size),
                    thrust::device_pointer_cast(b_ptr),
                    thrust::device_pointer_cast(res),
                    thrust::divides<float>());
}

void gpuCopy(float *v, float *dest, int size) {
  thrust::copy(thrust::device,
               thrust::device_pointer_cast(v),
               thrust::device_pointer_cast(v + size),
               thrust::device_pointer_cast(dest));
}

__device__
float det3x3(float *c, int id) {
  return c[0 + 9 * id] * c[4 + 9 * id] * c[8 + 9 * id] +
        c[1 + 9 * id] * c[6 + 9 * id] * c[5 + 9 * id] +
        c[2 + 9 * id] * c[3 + 9 * id] * c[7 + 9 * id] -
        c[6 + 9 * id] * c[4 + 9 * id] * c[2 + 9 * id] -
        c[0 + 9 * id] * c[5 + 9 * id] * c[7 + 9 * id] -
        c[1 + 9 * id] * c[3 + 9 * id] * c[8 + 9 * id];
}

__device__
float det33(float a0, float a1, float a2, float a3,
           float a4, float a5, float a6, float a7, float a8) {
  return a0 * a4 * a8 +
      a1 * a6 * a5 +
      a2 * a3 * a7 -
      a6 * a4 * a2 -
      a0 * a5 * a7 -
      a1 * a3 * a8;
}

__device__
float det4x4(float *c, int id) {
  float v1 = det33(c[5 + 16 * id], c[6 + 16 * id], c[7 + 16 * id], c[9 + 16 * id],
      c[10 + 16 * id], c[11 + 16 * id], c[13 + 16 * id], c[14 + 16 * id], c[15 + 16 * id]);
  float v2 = det33(c[1 + 16 * id], c[2 + 16 * id], c[3 + 16 * id], c[9 + 16 * id],
      c[10 + 16 * id], c[11 + 16 * id], c[13 + 16 * id], c[14 + 16 * id], c[15 + 16 * id]);
  float v3 = det33(c[1 + 16 * id], c[2 + 16 * id], c[3 + 16 * id], c[5 + 16 * id],
      c[6 + 16 * id], c[7 + 16 * id], c[13 + 16 * id], c[14 + 16 * id], c[15 + 16 * id]);
  float v4 = det33(c[1 + 16 * id], c[2 + 16 * id], c[3 + 16 * id], c[5 + 16 * id],
      c[6 + 16 * id], c[7 + 16 * id], c[9 + 16 * id], c[10 + 16 * id], c[11 + 16 * id]);
  return v1 - v2 + v3 - v4;
}

__device__
void inverse3x3(float *ic, float *c, float det, int id) {
  ic[0 + 3 * 0 + 9 * id] = (c[1 + 3 * 1 + 9 * id] * c[2 + 3 * 2 + 9 * id] - c[1 + 3 * 2 + 9 * id] * c[2 + 3 * 1 + 9 * id]) / det;
  ic[1 + 3 * 0 + 9 * id] = (c[2 + 3 * 0 + 9 * id] * c[1 + 3 * 2 + 9 * id] - c[1 + 3 * 0 + 9 * id] * c[2 + 3 * 2 + 9 * id]) / det;
  ic[2 + 3 * 0 + 9 * id] = (c[1 + 3 * 0 + 9 * id] * c[2 + 3 * 1 + 9 * id] - c[2 + 3 * 0 + 9 * id] * c[1 + 3 * 1 + 9 * id]) / det;

  ic[0 + 3 * 1 + 9 * id] = (c[2 + 3 * 1 + 9 * id] * c[0 + 3 * 2 + 9 * id] - c[0 + 3 * 1 + 9 * id] * c[2 + 3 * 2 + 9 * id]) / det;
  ic[1 + 3 * 1 + 9 * id] = (c[0 + 3 * 0 + 9 * id] * c[2 + 3 * 2 + 9 * id] - c[2 + 3 * 0 + 9 * id] * c[0 + 3 * 2 + 9 * id]) / det;
  ic[2 + 3 * 1 + 9 * id] = (c[0 + 3 * 1 + 9 * id] * c[2 + 3 * 0 + 9 * id] - c[0 + 3 * 0 + 9 * id] * c[2 + 3 * 1 + 9 * id]) / det;

  ic[0 + 3 * 2 + 9 * id] = (c[0 + 3 * 1 + 9 * id] * c[1 + 3 * 2 + 9 * id] - c[0 + 3 * 2 + 9 * id] * c[1 + 3 * 1 + 9 * id]) / det;
  ic[1 + 3 * 2 + 9 * id] = (c[0 + 3 * 2 + 9 * id] * c[1 + 3 * 0 + 9 * id] - c[0 + 3 * 0 + 9 * id] * c[1 + 3 * 2 + 9 * id]) / det;
  ic[2 + 3 * 2 + 9 * id] = (c[0 + 3 * 0 + 9 * id] * c[1 + 3 * 1 + 9 * id] - c[0 + 3 * 1 + 9 * id] * c[1 + 3 * 0 + 9 * id]) / det;
}

__device__
void inverse4x4(float *im, float *m, int id) {
  float inv[16], det;

  inv[0] = m[5 + 16 * id]  * m[10 + 16 * id] * m[15 + 16 * id] -
           m[5 + 16 * id]  * m[11 + 16 * id] * m[14 + 16 * id] -
           m[9 + 16 * id]  * m[6 + 16 * id]  * m[15 + 16 * id] +
           m[9 + 16 * id]  * m[7 + 16 * id]  * m[14 + 16 * id] +
           m[13 + 16 * id] * m[6 + 16 * id]  * m[11 + 16 * id] -
           m[13 + 16 * id] * m[7 + 16 * id]  * m[10 + 16 * id];

  inv[4] = -m[4 + 16 * id]  * m[10 + 16 * id] * m[15 + 16 * id] +
            m[4 + 16 * id]  * m[11 + 16 * id] * m[14 + 16 * id] +
            m[8 + 16 * id]  * m[6 + 16 * id]  * m[15 + 16 * id] -
            m[8 + 16 * id]  * m[7 + 16 * id]  * m[14 + 16 * id] -
            m[12 + 16 * id] * m[6 + 16 * id]  * m[11 + 16 * id] +
            m[12 + 16 * id] * m[7 + 16 * id]  * m[10 + 16 * id];

  inv[8] = m[4 + 16 * id]  * m[9 + 16 * id] * m[15 + 16 * id] -
           m[4 + 16 * id]  * m[11 + 16 * id] * m[13 + 16 * id] -
           m[8 + 16 * id]  * m[5 + 16 * id] * m[15 + 16 * id] +
           m[8 + 16 * id]  * m[7 + 16 * id] * m[13 + 16 * id] +
           m[12 + 16 * id] * m[5 + 16 * id] * m[11 + 16 * id] -
           m[12 + 16 * id] * m[7 + 16 * id] * m[9 + 16 * id];

  inv[12] = -m[4 + 16 * id]  * m[9 + 16 * id] * m[14 + 16 * id] +
             m[4 + 16 * id]  * m[10 + 16 * id] * m[13 + 16 * id] +
             m[8 + 16 * id]  * m[5 + 16 * id] * m[14 + 16 * id] -
             m[8 + 16 * id]  * m[6 + 16 * id] * m[13 + 16 * id] -
             m[12 + 16 * id] * m[5 + 16 * id] * m[10 + 16 * id] +
             m[12 + 16 * id] * m[6 + 16 * id] * m[9 + 16 * id];

  inv[1] = -m[1 + 16 * id]  * m[10 + 16 * id] * m[15 + 16 * id] +
            m[1 + 16 * id]  * m[11 + 16 * id] * m[14 + 16 * id] +
            m[9 + 16 * id]  * m[2 + 16 * id] * m[15 + 16 * id] -
            m[9 + 16 * id]  * m[3 + 16 * id] * m[14 + 16 * id] -
            m[13 + 16 * id] * m[2 + 16 * id] * m[11 + 16 * id] +
            m[13 + 16 * id] * m[3 + 16 * id] * m[10 + 16 * id];

  inv[5] = m[0 + 16 * id]  * m[10 + 16 * id] * m[15 + 16 * id] -
              m[0 + 16 * id]  * m[11 + 16 * id] * m[14 + 16 * id] -
              m[8 + 16 * id]  * m[2 + 16 * id] * m[15 + 16 * id] +
              m[8 + 16 * id]  * m[3 + 16 * id] * m[14 + 16 * id] +
              m[12 + 16 * id] * m[2 + 16 * id] * m[11 + 16 * id] -
              m[12 + 16 * id] * m[3 + 16 * id] * m[10 + 16 * id];

  inv[9] = -m[0 + 16 * id]  * m[9 + 16 * id] * m[15 + 16 * id] +
               m[0 + 16 * id]  * m[11 + 16 * id] * m[13 + 16 * id] +
               m[8 + 16 * id]  * m[1 + 16 * id] * m[15 + 16 * id] -
               m[8 + 16 * id]  * m[3 + 16 * id] * m[13 + 16 * id] -
               m[12 + 16 * id] * m[1 + 16 * id] * m[11 + 16 * id] +
               m[12 + 16 * id] * m[3 + 16 * id] * m[9 + 16 * id];

  inv[13] = m[0 + 16 * id]  * m[9 + 16 * id] * m[14 + 16 * id] -
               m[0 + 16 * id]  * m[10 + 16 * id] * m[13 + 16 * id] -
               m[8 + 16 * id]  * m[1 + 16 * id] * m[14 + 16 * id] +
               m[8 + 16 * id]  * m[2 + 16 * id] * m[13 + 16 * id] +
               m[12 + 16 * id] * m[1 + 16 * id] * m[10 + 16 * id] -
               m[12 + 16 * id] * m[2 + 16 * id] * m[9 + 16 * id];

  inv[2] = m[1 + 16 * id]  * m[6 + 16 * id] * m[15 + 16 * id] -
              m[1 + 16 * id]  * m[7 + 16 * id] * m[14 + 16 * id] -
              m[5 + 16 * id]  * m[2 + 16 * id] * m[15 + 16 * id] +
              m[5 + 16 * id]  * m[3 + 16 * id] * m[14 + 16 * id] +
              m[13 + 16 * id] * m[2 + 16 * id] * m[7 + 16 * id] -
              m[13 + 16 * id] * m[3 + 16 * id] * m[6 + 16 * id];

  inv[6] = -m[0 + 16 * id]  * m[6 + 16 * id] * m[15 + 16 * id] +
               m[0 + 16 * id]  * m[7 + 16 * id] * m[14 + 16 * id] +
               m[4 + 16 * id]  * m[2 + 16 * id] * m[15 + 16 * id] -
               m[4 + 16 * id]  * m[3 + 16 * id] * m[14 + 16 * id] -
               m[12 + 16 * id] * m[2 + 16 * id] * m[7 + 16 * id] +
               m[12 + 16 * id] * m[3 + 16 * id] * m[6 + 16 * id];

  inv[10] = m[0 + 16 * id]  * m[5 + 16 * id] * m[15 + 16 * id] -
               m[0 + 16 * id]  * m[7 + 16 * id] * m[13 + 16 * id] -
               m[4 + 16 * id]  * m[1 + 16 * id] * m[15 + 16 * id] +
               m[4 + 16 * id]  * m[3 + 16 * id] * m[13 + 16 * id] +
               m[12 + 16 * id] * m[1 + 16 * id] * m[7 + 16 * id] -
               m[12 + 16 * id] * m[3 + 16 * id] * m[5 + 16 * id];

  inv[14] = -m[0 + 16 * id]  * m[5 + 16 * id] * m[14 + 16 * id] +
                m[0 + 16 * id]  * m[6 + 16 * id] * m[13 + 16 * id] +
                m[4 + 16 * id]  * m[1 + 16 * id] * m[14 + 16 * id] -
                m[4 + 16 * id]  * m[2 + 16 * id] * m[13 + 16 * id] -
                m[12 + 16 * id] * m[1 + 16 * id] * m[6 + 16 * id] +
                m[12 + 16 * id] * m[2 + 16 * id] * m[5 + 16 * id];

  inv[3] = -m[1 + 16 * id] * m[6 + 16 * id] * m[11 + 16 * id] +
               m[1 + 16 * id] * m[7 + 16 * id] * m[10 + 16 * id] +
               m[5 + 16 * id] * m[2 + 16 * id] * m[11 + 16 * id] -
               m[5 + 16 * id] * m[3 + 16 * id] * m[10 + 16 * id] -
               m[9 + 16 * id] * m[2 + 16 * id] * m[7 + 16 * id] +
               m[9 + 16 * id] * m[3 + 16 * id] * m[6 + 16 * id];

  inv[7] = m[0 + 16 * id] * m[6 + 16 * id] * m[11 + 16 * id] -
              m[0 + 16 * id] * m[7 + 16 * id] * m[10 + 16 * id] -
              m[4 + 16 * id] * m[2 + 16 * id] * m[11 + 16 * id] +
              m[4 + 16 * id] * m[3 + 16 * id] * m[10 + 16 * id] +
              m[8 + 16 * id] * m[2 + 16 * id] * m[7 + 16 * id] -
              m[8 + 16 * id] * m[3 + 16 * id] * m[6 + 16 * id];

  inv[11] = -m[0 + 16 * id] * m[5 + 16 * id] * m[11 + 16 * id] +
                m[0 + 16 * id] * m[7 + 16 * id] * m[9 + 16 * id] +
                m[4 + 16 * id] * m[1 + 16 * id] * m[11 + 16 * id] -
                m[4 + 16 * id] * m[3 + 16 * id] * m[9 + 16 * id] -
                m[8 + 16 * id] * m[1 + 16 * id] * m[7 + 16 * id] +
                m[8 + 16 * id] * m[3 + 16 * id] * m[5 + 16 * id];

  inv[15] = m[0 + 16 * id] * m[5 + 16 * id] * m[10 + 16 * id] -
               m[0 + 16 * id] * m[6 + 16 * id] * m[9 + 16 * id] -
               m[4 + 16 * id] * m[1 + 16 * id] * m[10 + 16 * id] +
               m[4 + 16 * id] * m[2 + 16 * id] * m[9 + 16 * id] +
               m[8 + 16 * id] * m[1 + 16 * id] * m[6 + 16 * id] -
               m[8 + 16 * id] * m[2 + 16 * id] * m[5 + 16 * id];

     det = m[0 + 16 * id] * inv[0] + m[1 + 16 * id] * inv[4] + m[2 + 16 * id] * inv[8] + m[3 + 16 * id] * inv[12];
     if (det == 0.0f)
       printf("DET == 0!!! -> %d\n", id);
     det = 1.0 / det;

     for (int i = 0; i < 16; i++)
         im[i + 16 * id] = inv[i] * det;
}

__device__
float det(float *c, int size) {
  if (size == 1) {
    return c[0];
  } else if (size == 2) {
    return c[0 + 0 * 2] * c[1 + 1 * 2] - c[0 + 1 * 2] * c[1 + 0 * 2];
  } else if (size == 3) {
    return c[0]*c[4]*c[8] +
        c[1]*c[6]*c[5] +
        c[2]*c[3]*c[7] -
        c[6]*c[4]*c[2] -
        c[0]*c[5]*c[7] -
        c[1]*c[3]*c[8];
  } else if (size == 4) {
    float v1 = det33(c[5], c[6], c[7], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v2 = det33(c[1], c[2], c[3], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v3 = det33(c[1], c[2], c[3], c[5], c[6], c[7], c[13], c[14], c[15]);
    float v4 = det33(c[1], c[2], c[3], c[5], c[6], c[7], c[9], c[10], c[11]);
    return v1 - v2 + v3 - v4;
  }
}

void gpuReductionWithMaskAndTransform(float *v, float *mask, int size, float *res, int size_new) {
  CheckRunTime(__func__)
  thrust::device_vector<float> d_v(v, v + size);
  thrust::device_vector<float> d_mask(mask, mask + size);

  thrust::device_vector<float> d_res(size_new);
  thrust::device_vector<float> d_mask_new(size_new);

  thrust::sort_by_key(d_mask.begin(), d_mask.end(), d_v.begin());
  thrust::reduce_by_key(d_mask.begin(), d_mask.end(), d_v.begin(), d_mask_new.begin(), d_res.begin());

  for (int i = 0; i < size; ++i) {
    float tmp = d_res[static_cast<int>(mask[i])];
    res[i] = tmp;
    //std::cout << tmp << " ";
  }
  //    exit(-1);
}

__global__ void kernelTransformFrom_2N_to_6E(float *mask, float *a, float *b, int size) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    b[i] = a[(int)mask[i]];
  }
}

void gpuTransformFrom_2N_to_6E(float *mask_ptr, float *a_ptr,
                               float *b_ptr, int size) {
  //    float* mask_ptr = thrust::raw_pointer_cast( mask_vec.data() );
  //    float* a_ptr = thrust::raw_pointer_cast( a_vec.data() );
  //    float* b_ptr = thrust::raw_pointer_cast( b_vec.data() );
  kernelTransformFrom_2N_to_6E<<<(size + 255) / 256, 256>>>(mask_ptr, a_ptr, b_ptr, size);
}

void gpuReductionWithMaskAndTransform2(float *v, float *mask, int size,
                                       float *res, int size_new) {
  CheckRunTime(__func__)
  thrust::device_vector<float> temp_res(size_new);
  thrust::device_vector<float> temp_mask(size);
  thrust::device_vector<float> temp_v(size);

  thrust::copy(thrust::device, mask, mask + size, temp_mask.begin());
  thrust::copy(thrust::device, v, v + size, temp_v.begin());

  thrust::sort_by_key(temp_mask.begin(), temp_mask.end(), temp_v.begin());
  thrust::reduce_by_key(temp_mask.begin(), temp_mask.end(), temp_v.begin(), thrust::make_discard_iterator(), temp_res.begin());

  gpuTransformFrom_2N_to_6E(mask, thrust::raw_pointer_cast(temp_res.data()), res, size);
}

void gpuReductionWithMask(float *v, float *mask, int size, float *res, int size_new) {
  CheckRunTime(__func__)
  thrust::device_vector<float> d_v(v, v + size);
  thrust::device_vector<float> d_mask(mask, mask + size);

  thrust::device_vector<float> d_res(size_new);

  thrust::sort_by_key(d_mask.begin(), d_mask.end(), d_v.begin());
  thrust::reduce_by_key(d_mask.begin(), d_mask.end(), d_v.begin(), thrust::make_discard_iterator(), d_res.begin());

  for (int i = 0; i < size_new; ++i) {
    float tmp = d_res[i];
    res[i] = tmp;
  }
}

void gpuReductionWithMask2(float *v, float *mask,
                           int size, float *res) {
  CheckRunTime(__func__)
  thrust::device_vector<float> temp_v(size);
  thrust::device_vector<float> temp_mask(size);
  thrust::copy(thrust::device, thrust::device_pointer_cast(v),
               thrust::device_pointer_cast(v + size), temp_v.begin());

  thrust::copy(thrust::device, thrust::device_pointer_cast(mask),
               thrust::device_pointer_cast(mask + size), temp_mask.begin());
  thrust::sort_by_key(temp_mask.begin(), temp_mask.end(), temp_v.begin());
  thrust::reduce_by_key(temp_mask.begin(), temp_mask.end(), temp_v.begin(),
                        thrust::make_discard_iterator(), thrust::device_pointer_cast(res));
}

void thrustReductionWithMask(thrust::device_vector<float> &d_v,
                             thrust::device_vector<float> &d_mask,
                             thrust::device_vector<float> &d_res) {
    CheckRunTime(__func__)
    thrust::device_vector<float> mask_sorted(d_mask.size()), v_sorted(d_v.size());

    mask_sorted = d_mask;
    v_sorted = d_v;
    thrust::sort_by_key(mask_sorted.begin(), mask_sorted.end(), v_sorted.begin());

    thrust::reduce_by_key(mask_sorted.begin(), mask_sorted.end(), v_sorted.begin(),
                          thrust::make_discard_iterator(), d_res.begin());
}

void thrustCountNAdjElem(thrust::device_vector<float> &d_mask,
                         thrust::device_vector<float> &d_res) {
    CheckRunTime(__func__)
    thrust::device_vector<float> mask_sorted(d_mask.size());

    mask_sorted = d_mask;
    thrust::sort(mask_sorted.begin(), mask_sorted.end());

    thrust::reduce_by_key(mask_sorted.begin(), mask_sorted.end(), thrust::make_constant_iterator(1.f),
                          thrust::make_discard_iterator(), d_res.begin());
}

void gpuTransform_2N_to_6E(float *d_v, int n_gl_dofs, float *d_mask, float *d_res, int grid_size) {
    CheckRunTime(__func__)
    // Make parallel!
    thrust::host_vector<float> host_res(grid_size), host_v(n_gl_dofs), host_mask(grid_size);
    gpuCopyDeviceToHost(d_v, thrust::raw_pointer_cast(host_v.data()), n_gl_dofs);
    //gpuCopy(d_v, thrust::raw_pointer_cast(host_v.data()), n_gl_dofs);
    gpuCopyDeviceToHost(d_mask, thrust::raw_pointer_cast(host_mask.data()), grid_size);
    //gpuCopy(d_mask, thrust::raw_pointer_cast(host_mask.data()), grid_size);
    for (int i = 0; i < grid_size; ++i) {
      //TODO: CUDA
        host_res[i] = host_v[static_cast<int>(host_mask[i])];
    }

    gpuCopyHostToDevice(thrust::raw_pointer_cast(host_res.data()), d_res, grid_size);
    //gpuCopy(thrust::raw_pointer_cast(host_res.data()), d_res, grid_size);
//    cudaMemcpy(d_res, thrust::raw_pointer_cast(host_res.data()), grid_size * sizeof(float), cudaMemcpyHostToDevice);

}

void thrustTransform_2N_to_6E(thrust::device_vector<float> &d_v,
                              thrust::device_vector<float> &d_mask,
                              thrust::device_vector<float> &d_res) {
    CheckRunTime(__func__)
    // Make parallel!
    thrust::host_vector<float> host_res(d_res.size()), host_v(d_v.size()), host_mask(d_mask.size());
    host_v = d_v;
    host_mask = d_mask;
    for (int i = 0; i < d_mask.size(); ++i) {
        host_res[i] = host_v[static_cast<int>(host_mask[i])];
    }

    d_res = host_res;

}

float gpuDotProduct(float *a, float *b, int size) {
  CheckRunTime(__func__)
  thrust::device_vector<float> d_a(a, a + size);
  thrust::device_vector<float> d_b(b, b + size);
  return thrust::inner_product(thrust::cuda::par, d_a.begin(), d_a.end(), d_b.begin(), 0.0f);
}

float gpuDotProduct2(float *a, float *b, int size) {
  CheckRunTime(__func__)
  return thrust::inner_product(thrust::cuda::par, thrust::device_pointer_cast(a),
                               thrust::device_pointer_cast(a + size), thrust::device_pointer_cast(b), 0.0f);
}

void gpuDivideByElementwise(float *a, float *b,
                            float *c, int size) {
  CheckRunTime(__func__)
//      thrust::transform(a.begin(), a.end(), b.begin(), c.begin(),
//                        thrust::placeholders::_1 / thrust::placeholders::_2);
  thrust::transform(thrust::cuda::par, thrust::device_pointer_cast(a),
                    thrust::device_pointer_cast(a + size),
                    thrust::device_pointer_cast(b),
                    thrust::device_pointer_cast(c),
                    thrust::divides<float>());
}

__device__ void matrixMultiply(float *a, int row, int n, int col, float *b, float *res, int id) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      res[j + i * n + row * col * id] = 0.0f;
      for (int k = 0; k < n; k++) {
        res[j + i * n + row * col * id] += a[k + i * n + row * n * id] * b[j + k * col];
      }
    }
  }
}

__global__ void kernelCalculateKlocal2D(int elementsCount, float *elements, float *diag,
                                       float *nodesX, float *nodesY,
                                       float *D, float *K, float *B,
                                       float *constraints, int constraintsCount,
                                       float *C, float *IC, float *temp, float *B_T) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {

    C[0 + 0 * 3 + 9 * index] = C[0 + 1 * 3 + 9 * index] = C[0 + 2 * 3 + 9 * index] = 1.0f;
    C[1 + 0 * 3 + 9 * index] = nodesX[int(elements[3 * index + 0])];
    C[1 + 1 * 3 + 9 * index] = nodesX[int(elements[3 * index + 1])];
    C[1 + 2 * 3 + 9 * index] = nodesX[int(elements[3 * index + 2])];
    C[2 + 0 * 3 + 9 * index] = nodesY[int(elements[3 * index + 0])];
    C[2 + 1 * 3 + 9 * index] = nodesY[int(elements[3 * index + 1])];
    C[2 + 2 * 3 + 9 * index] = nodesY[int(elements[3 * index + 2])];

    float determinant = det3x3(C, index);
    inverse3x3(IC, C, determinant, index);

    for (int i = 0; i < 3; i++) {
      B[(2 * i + 0) + 0 * 6 + index * 18] = IC[i + 3 * 1 + 9 * index];
      B[(2 * i + 1) + 0 * 6 + index * 18] = 0.0f;
      B[(2 * i + 0) + 1 * 6 + index * 18] = 0.0f;
      B[(2 * i + 1) + 1 * 6 + index * 18] = IC[i + 3 * 2 + 9 * index];
      B[(2 * i + 0) + 2 * 6 + index * 18] = IC[i + 3 * 2 + 9 * index];
      B[(2 * i + 1) + 2 * 6 + index * 18] = IC[i + 3 * 1 + 9 * index];
    }

    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++) {
        B_T[j + i * 3 + 18 * index] = B[i + j * 6 + index * 18];
      }
    }

    //temp(6,3) = B_T(6,3) * D(3,3)
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++) {
        temp[j + i * 3 + 18 * index] = 0.0f;
        for (int k = 0; k < 3; k++) {
          temp[j + i * 3 + 18 * index] += B_T[k + i * 3 + 18 * index] * D[j + k * 3];
        }
//        if (index == 1)
//          printf("%f ", temp[j + i * 3 + 18 * index]);
      }
    }
//    matrixMultiply(B_T, 6, 3, 3, D, temp, index);

    //K(6,6) = temp(6,3) * B(3,6)
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        K[j + i * 6 + 36 * index] = 0.0f;
        for (int k = 0; k < 3; k++) {
          K[j + i * 6 + 36 * index] += temp[k + i * 3 + 6 * 3 * index] * B[j + k * 6 + index * 3 * 6];
        }
      }
    }
    //matrixMultiply(temp, 6, 3, 6, B, K, index);

    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        K[j + i * 6 + 36 * index] *= fabs(determinant) * 0.5f;
      }
    }

    for (int c_id = 0; c_id < constraintsCount; ++c_id) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int ilocal = 0; ilocal < 2; ++ilocal) {
            for (int jlocal = 0; jlocal < 2; ++jlocal) {
              if (2 * int(elements[3 * index + i]) + ilocal == int(constraints[c_id]) ||
                  2 * int(elements[3 * index + j]) + jlocal == int(constraints[c_id])) {
                if (2 * int(elements[3 * index + i]) + ilocal != 2 * int(elements[3 * index + j]) + jlocal) {
                  K[(2 * j + jlocal) + 6 * (2 * i + ilocal) + 36 * index] = 0.0f;
                }
              }
            }
          }
        }
      }
    }

    // Grisha: Calculate separately!
    diag[6 * index + 0] = K[0 + 6 * 0 + 36 * index];
    diag[6 * index + 1] = K[1 + 6 * 1 + 36 * index];
    diag[6 * index + 2] = K[2 + 6 * 2 + 36 * index];
    diag[6 * index + 3] = K[3 + 6 * 3 + 36 * index];
    diag[6 * index + 4] = K[4 + 6 * 4 + 36 * index];
    diag[6 * index + 5] = K[5 + 6 * 5 + 36 * index];
  }

  __syncthreads();

}

__global__ void kernelCalculateKlocal3D(int elementsCount, float *elements, float *diag,
                                       float *nodesX, float *nodesY, float *nodesZ,
                                       float *D, float *K, float *B,
                                       float *constraints, int constraintsCount,
                                       float *C, float *IC, float *temp, float *B_T) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {

    C[0 + 0 * 4 + 16 * index] = C[0 + 1 * 4 + 16 * index] = C[0 + 2 * 4 + 16 * index] = C[0 + 3 * 4 + 16 * index] = 1.0f;

    C[1 + 0 * 4 + 16 * index] = nodesX[int(elements[4 * index + 0])];
    C[1 + 1 * 4 + 16 * index] = nodesX[int(elements[4 * index + 1])];
    C[1 + 2 * 4 + 16 * index] = nodesX[int(elements[4 * index + 2])];
    C[1 + 3 * 4 + 16 * index] = nodesX[int(elements[4 * index + 3])];

    C[2 + 0 * 4 + 16 * index] = nodesY[int(elements[4 * index + 0])];
    C[2 + 1 * 4 + 16 * index] = nodesY[int(elements[4 * index + 1])];
    C[2 + 2 * 4 + 16 * index] = nodesY[int(elements[4 * index + 2])];
    C[2 + 3 * 4 + 16 * index] = nodesY[int(elements[4 * index + 3])];

    C[3 + 0 * 4 + 16 * index] = nodesZ[int(elements[4 * index + 0])];
    C[3 + 1 * 4 + 16 * index] = nodesZ[int(elements[4 * index + 1])];
    C[3 + 2 * 4 + 16 * index] = nodesZ[int(elements[4 * index + 2])];
    C[3 + 3 * 4 + 16 * index] = nodesZ[int(elements[4 * index + 3])];

    float determinant = det4x4(C, index);
    inverse4x4(IC, C, index);

    if (0) {
      printf("DET = %f\n", determinant);

      for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            printf("%f ", C[j + i * 4]);
          }
          printf("\n");
        }
        printf("\n");

        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            printf("%f ", IC[j + i * 4]);
          }
          printf("\n");
        }

//      printf("\n");
//      for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//          float temp = 0.0f;
//          for (int k = 0; k < 4; k++) {
//            temp += C[k + i * 4] * IC[j + k * 4];
//          }
//          printf("%f ", temp);
//        }
//        printf("\n");
//      }
    }

    for (int i = 0; i < 4; i++) {
      B[(3 * i + 0) + 0 * 12 + index * 6 * 12] = IC[i + 1 * 4 + 16 * index];
      B[(3 * i + 1) + 0 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 2) + 0 * 12 + index * 6 * 12] = 0.0f;

      B[(3 * i + 0) + 1 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 1) + 1 * 12 + index * 6 * 12] = IC[i + 2 * 4 + 16 * index];
      B[(3 * i + 2) + 1 * 12 + index * 6 * 12] = 0.0f;

      B[(3 * i + 0) + 2 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 1) + 2 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 2) + 2 * 12 + index * 6 * 12] = IC[i + 3 * 4 + 16 * index];

      B[(3 * i + 0) + 3 * 12 + index * 6 * 12] = IC[i + 2 * 4 + 16 * index];
      B[(3 * i + 1) + 3 * 12 + index * 6 * 12] = IC[i + 1 * 4 + 16 * index];
      B[(3 * i + 2) + 3 * 12 + index * 6 * 12] = 0.0f;

      B[(3 * i + 0) + 4 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 1) + 4 * 12 + index * 6 * 12] = IC[i + 3 * 4 + 16 * index];
      B[(3 * i + 2) + 4 * 12 + index * 6 * 12] = IC[i + 2 * 4 + 16 * index];

      B[(3 * i + 0) + 5 * 12 + index * 6 * 12] = IC[i + 3 * 4 + 16 * index];
      B[(3 * i + 1) + 5 * 12 + index * 6 * 12] = 0.0f;
      B[(3 * i + 2) + 5 * 12 + index * 6 * 12] = IC[i + 1 * 4 + 16 * index];
    }

    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 6; j++) {
        B_T[j + i * 6 + 6 * 12 * index] = B[i + j * 12 + index * 6 * 12];
      }
    }

    //temp(12,6) = B_T(12,6) * D(6,6)
    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 6; j++) {
        temp[j + i * 6 + 6 * 12 * index] = 0.0f;
        for (int k = 0; k < 6; k++) {
          temp[j + i * 6 + 6 * 12 * index] += B_T[k + i * 6 + 6 * 12 * index] * D[j + k * 6];
        }
      }
    }
//    matrixMultiply(B_T, 6, 3, 3, D, temp, index);

    //K(12,12) = temp(12,6) * B(6,12)
    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 12; j++) {
        K[j + i * 12 + 12 * 12 * index] = 0.0f;
        for (int k = 0; k < 6; k++) {
          K[j + i * 12 + 12 * 12 * index] += temp[k + i * 6 + 12 * 6 * index] * B[j + k * 12 + index * 6 * 12];
        }
      }
    }
    //matrixMultiply(temp, 6, 3, 6, B, K, index);

    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 12; j++) {
        K[j + i * 12 + 12 * 12 * index] *= fabs(determinant) / 6.f;
      }
    }

    for (int c_id = 0; c_id < constraintsCount; ++c_id) {
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          for (int ilocal = 0; ilocal < 3; ++ilocal) {
            for (int jlocal = 0; jlocal < 3; ++jlocal) {
              if (3 * int(elements[4 * index + i]) + ilocal == int(constraints[c_id]) ||
                  3 * int(elements[4 * index + j]) + jlocal == int(constraints[c_id])) {
                if (3 * int(elements[4 * index + i]) + ilocal != 3 * int(elements[4 * index + j]) + jlocal) {
                  K[(3 * j + jlocal) + 12 * (3 * i + ilocal) + 12 * 12 * index] = 0.0f;
                }
              }
            }
          }
        }
      }
    }

    // Grisha: Calculate separately!
    diag[12 * index + 0] = K[0 + 12 * 0 + 12 * 12 * index];
    diag[12 * index + 1] = K[1 + 12 * 1 + 12 * 12 * index];
    diag[12 * index + 2] = K[2 + 12 * 2 + 12 * 12 * index];
    diag[12 * index + 3] = K[3 + 12 * 3 + 12 * 12 * index];
    diag[12 * index + 4] = K[4 + 12 * 4 + 12 * 12 * index];
    diag[12 * index + 5] = K[5 + 12 * 5 + 12 * 12 * index];
    diag[12 * index + 6] = K[6 + 12 * 6 + 12 * 12 * index];
    diag[12 * index + 7] = K[7 + 12 * 7 + 12 * 12 * index];
    diag[12 * index + 8] = K[8 + 12 * 8 + 12 * 12 * index];
    diag[12 * index + 9] = K[9 + 12 * 9 + 12 * 12 * index];
    diag[12 * index + 10] = K[10 + 12 * 10 + 12 * 12 * index];
    diag[12 * index + 11] = K[11 + 12 * 11 + 12 * 12 * index];
  }

  __syncthreads();

}

void gpuCalculateKlocal2(gpuDataKeeper &gpu_data, FEMdataKeeper &FEMdata) {
  CheckRunTime(__func__)

  int DIM = FEMdata.DIM;

  int elementsCount = FEMdata.elementsCount;
  int nodesCount = FEMdata.nodesCount;
  int constraintsCount = FEMdata.CudaIndicesToConstraintsCount;

  float *h_D = FEMdata.D.get_data();
  float *h_constraints = FEMdata.CudaIndicesToConstraints.get_data();
  float *h_nodesX = FEMdata.nodes[0].get_data();
  float *h_nodesY = FEMdata.nodes[1].get_data();

  float* elements = gpu_data.get_Elements();//thrust::raw_pointer_cast( CudaElements.data() );
  float* diag = gpu_data.get_diag();//thrust::raw_pointer_cast( diag_vec.data() );

  float *d_constraints;
  float *d_nodesX, *d_nodesY, *d_D;

  cudaMalloc((void**)&d_D, 3 * (DIM - 1) * 3 * (DIM - 1) * sizeof(float));
  cudaMemcpy(d_D, h_D, 3 * (DIM - 1) * 3 * (DIM - 1) * sizeof(float), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_constraints, constraintsCount * sizeof(float));
  cudaMemcpy(d_constraints, h_constraints, constraintsCount * sizeof(float), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_nodesX, nodesCount * sizeof(float));
  cudaMalloc((void**)&d_nodesY, nodesCount * sizeof(float));
  cudaMemcpy(d_nodesX, h_nodesX, nodesCount * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nodesY, h_nodesY, nodesCount * sizeof(float), cudaMemcpyHostToDevice);

  float *d_C, *d_IC, *d_temp, *d_B_T;
  cudaMalloc((void**)&d_C, (DIM + 1) * (DIM + 1) * elementsCount * sizeof(float));
  cudaMalloc((void**)&d_IC, (DIM + 1) * (DIM + 1) * elementsCount * sizeof(float));
  cudaMalloc((void**)&d_temp, 3 * (DIM - 1) * 6 * (DIM - 1) * elementsCount * sizeof(float));
  cudaMalloc((void**)&d_B_T, 3 * (DIM - 1) * 6 * (DIM - 1) * elementsCount * sizeof(float));

  if (DIM == 2) {
    kernelCalculateKlocal2D<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements, diag,
                                                               d_nodesX, d_nodesY,
                                                               d_D, gpu_data.get_Klocals(), gpu_data.get_B(),
                                                               d_constraints, constraintsCount, d_C, d_IC, d_temp, d_B_T);
  } else if (DIM == 3) {
    float *h_nodesZ = FEMdata.nodes[2].get_data();
    float *d_nodesZ;
    cudaMalloc((void**)&d_nodesZ, nodesCount * sizeof(float));
    cudaMemcpy(d_nodesZ, h_nodesZ, nodesCount * sizeof(float), cudaMemcpyHostToDevice);

    kernelCalculateKlocal3D<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements, diag,
                                                               d_nodesX, d_nodesY, d_nodesZ,
                                                               d_D, gpu_data.get_Klocals(), gpu_data.get_B(),
                                                               d_constraints, constraintsCount, d_C, d_IC, d_temp, d_B_T);
    cudaFree(d_nodesZ);
  }


  cudaFree(d_nodesX);
  cudaFree(d_nodesY);
  cudaFree(d_constraints);

  cudaFree(d_C);
  cudaFree(d_IC);
  cudaFree(d_temp);
  cudaFree(d_B_T);
}

__global__ void kernelCalculateMlocal(int elementsCount, float *elements,
                                       float *nodesX, float *nodesY, float *M, float rho, bool isLumped) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    float X0 = nodesX[int(elements[3 * index + 0])], X1 = nodesX[int(elements[3 * index + 1])], X2 = nodesX[int(elements[3 * index + 2])];
    float Y0 = nodesY[int(elements[3 * index + 0])], Y1 = nodesY[int(elements[3 * index + 1])], Y2 = nodesY[int(elements[3 * index + 2])];
    float area = 0.5 * std::abs( (X0 - X2) * (Y1 - Y2) - (X1 - X2) * (Y0 - Y2) );
    float mass = rho * area;

    if (isLumped) {
      M[0 + 6 * index] = mass / 3;
      M[1 + 6 * index] = mass / 3;
      M[2 + 6 * index] = mass / 3;
      M[3 + 6 * index] = mass / 3;
      M[4 + 6 * index] = mass / 3;
      M[5 + 6 * index] = mass / 3;
//      M[0 + 0 * 6 + 36 * index] = mass / 3;
//      M[1 + 1 * 6 + 36 * index] = mass / 3;
//      M[2 + 2 * 6 + 36 * index] = mass / 3;
//      M[3 + 3 * 6 + 36 * index] = mass / 3;
//      M[4 + 4 * 6 + 36 * index] = mass / 3;
//      M[5 + 5 * 6 + 36 * index] = mass / 3;
    } else {
      M[0 + 0 * 6 + 36 * index] = mass / 6;
      M[1 + 1 * 6 + 36 * index] = mass / 6;
      M[2 + 2 * 6 + 36 * index] = mass / 6;
      M[3 + 3 * 6 + 36 * index] = mass / 6;
      M[4 + 4 * 6 + 36 * index] = mass / 6;
      M[5 + 5 * 6 + 36 * index] = mass / 6;

      M[0 + 2 * 6 + 36 * index] = mass / 12;
      M[1 + 3 * 6 + 36 * index] = mass / 12;
      M[0 + 4 * 6 + 36 * index] = mass / 12;
      M[1 + 5 * 6 + 36 * index] = mass / 12;
      M[2 + 4 * 6 + 36 * index] = mass / 12;
      M[3 + 5 * 6 + 36 * index] = mass / 12;
      M[2 + 0 * 6 + 36 * index] = mass / 12;
      M[3 + 1 * 6 + 36 * index] = mass / 12;
      M[4 + 0 * 6 + 36 * index] = mass / 12;
      M[5 + 1 * 6 + 36 * index] = mass / 12;
      M[4 + 2 * 6 + 36 * index] = mass / 12;
      M[5 + 3 * 6 + 36 * index] = mass / 12;
    }

  }

}

// 2D only
void gpuCalculateMlocal(gpuDataKeeper_DYN &gpu_data, FEMdataKeeper &FEMdata, float rho) {
  CheckRunTime(__func__)
  float* elements = gpu_data.get_Elements();

  float *d_nodesX, *d_nodesY;
  cudaMalloc((void**)&d_nodesX, FEMdata.nodesCount * sizeof(float));
  cudaMalloc((void**)&d_nodesY, FEMdata.nodesCount * sizeof(float));
  cudaMemcpy(d_nodesX, FEMdata.nodes[0].get_data(), FEMdata.nodesCount * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nodesY, FEMdata.nodes[1].get_data(), FEMdata.nodesCount * sizeof(float), cudaMemcpyHostToDevice);

  if (gpu_data.isLumped) {
    kernelCalculateMlocal<<<(FEMdata.elementsCount + 255) / 256, 256>>>(FEMdata.elementsCount, elements,
                                                                d_nodesX, d_nodesY, gpu_data.get_diagM(), rho, true);
//    kernelCalculateMlocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements,
//                                                                d_nodesX, d_nodesY, gpu_data.get_Mlocals(), rho, true);
  } else {
    kernelCalculateMlocal<<<(FEMdata.elementsCount + 255) / 256, 256>>>(FEMdata.elementsCount, elements,
                                                                d_nodesX, d_nodesY, gpu_data.get_Mlocals(), rho, false);
  }
  cudaFree(d_nodesX);
  cudaFree(d_nodesY);
}

__global__ void kernelMultiply(int elementsCount, float *u, float *Klocals, float *p) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    for (int j = 0; j < 6; j++) {
      u[j + index * 6] = 0.0f;
      for (int k = 0; k < 6; k++) {
        u[j + index * 6] += Klocals[k + j * 6 + index * 36] * p[k + index * 6];
      }
    }
  }
}

__global__ void kernelMultiply3D(int elementsCount, float *u, float *Klocals, float *p) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    for (int j = 0; j < 12; j++) {
      u[j + index * 12] = 0.0f;
      for (int k = 0; k < 12; k++) {
        u[j + index * 12] += Klocals[k + j * 12 + index * 12 * 12] * p[k + index * 12];
      }
    }
  }
}

__global__ void kernelMultiply_2WeightedMatrix(int elementsCount, float *u, float *Mlocals, float *Klocals, float *p,
                                               float cM, float cK) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    for (int j = 0; j < 6; j++) {
      u[j + index * 6] = 0.0f;
      for (int k = 0; k < 6; k++) {
        u[j + index * 6] +=
            (cM*Mlocals[k + j * 6 + index * 36] + cK*Klocals[k + j * 6 + index * 36]) *
            p[k + index * 6];
      }
    }
  }
}

__global__ void kernelMultiply_DAMP(int elementsCount, float *u, float *Mlocals, float *Klocals, float *p,
                                    float cM, float cK, float cC, float damping_alpha, float damping_beta, float isLumped) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    if (isLumped) {
      for (int j = 0; j < 6; j++) {
        u[j + index * 6] = 0.0f;
        for (int k = 0; k < 6; k++) {
          float factor = (cK + cC*damping_beta)*Klocals[k + j * 6 + index * 36];
          if (k == j) factor += (cM + cC*damping_alpha)*Mlocals[k + index * 6];
          u[j + index * 6] += factor * p[k + index * 6];
        }
      }
    } else {
      for (int j = 0; j < 6; j++) {
        u[j + index * 6] = 0.0f;
        for (int k = 0; k < 6; k++) {
          u[j + index * 6] +=
              (cM*Mlocals[k + j * 6 + index * 36] + cK*Klocals[k + j * 6 + index * 36] +
              cC*(damping_alpha*Mlocals[k + j * 6 + index * 36] + damping_beta*Klocals[k + j * 6 + index * 36])) *
              p[k + index * 6];
       }
      }
    }
  }
}

__global__ void kernelMultiply_Clocal(int elementsCount, float *u, float *Mlocals, float *Klocals, float *p,
                                      float damping_alpha, float damping_beta, bool isLumped) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    if (isLumped) {
      for (int j = 0; j < 6; j++) {
        u[j + index * 6] = 0.0f;
        for (int k = 0; k < 6; k++) {
          // if lumped (explicit scheme), then damping_beta = 0 !!!
          //float factor = 0.0f;
          float factor = damping_beta*Klocals[k + j * 6 + index * 36];
          if (k == j)
            factor += damping_alpha*Mlocals[k + index * 6];
          u[j + index * 6] += factor * p[k + index * 6];
        }
      }
    } else {
      for (int j = 0; j < 6; j++) {
        u[j + index * 6] = 0.0f;
        for (int k = 0; k < 6; k++) {
          u[j + index * 6] +=
              (damping_alpha*Mlocals[k + j * 6 + index * 36] + damping_beta*Klocals[k + j * 6 + index * 36]) *
              p[k + index * 6];
        }
    }
    }
  }
}


void gpuMultiplyKlocalByVec(gpuDataKeeper &gpu_data, int DIM, int elementsCount) {
  CheckRunTime(__func__)
  float* u_ptr = gpu_data.get_u();//thrust::raw_pointer_cast( u.data() );
  float* p_ptr = gpu_data.get_p();//thrust::raw_pointer_cast( p.data() );
  if (DIM == 2) {
    kernelMultiply<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Klocals(), p_ptr);
  } else if (DIM == 3) {
    kernelMultiply3D<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Klocals(), p_ptr);
  }
}

void gpuMultiplyMatrixByVec(float* Matr, float* Vec, float* Res, int elementsCount) {
  CheckRunTime(__func__)
  kernelMultiply<<<(elementsCount + 255) / 256, 256>>>(elementsCount, Res, Matr, Vec);
}

void gpuMultiplyClocalByVec(gpuDataKeeper_DYN_DAMP &gpu_data, float* Vec, float* Res, int elementsCount) {
  CheckRunTime(__func__)
  if (gpu_data.isLumped) {
    kernelMultiply_Clocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, Res, gpu_data.get_diagM(), gpu_data.get_Klocals(), Vec,
                                                                gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, true);
  } else {
    kernelMultiply_Clocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, Res, gpu_data.get_Mlocals(), gpu_data.get_Klocals(), Vec,
                                                                gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, false);
  }
}

void gpuMultiplyAlocalByVec(gpuDataKeeper_DYN &gpu_data, int elementsCount) {
  CheckRunTime(__func__)
  float* u_ptr = gpu_data.get_u();//thrust::raw_pointer_cast( u.data() );
  float* p_ptr = gpu_data.get_p();//thrust::raw_pointer_cast( p.data() );
  kernelMultiply_2WeightedMatrix<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Mlocals(), gpu_data.get_Klocals(), p_ptr,
                                                                       gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2);
}

void gpuMultiplyAlocalByVec_DAMP(gpuDataKeeper_DYN_DAMP &gpu_data, int elementsCount) {
  CheckRunTime(__func__)
  float* u_ptr = gpu_data.get_u();//thrust::raw_pointer_cast( u.data() );
  float* p_ptr = gpu_data.get_p();//thrust::raw_pointer_cast( p.data() );
  if (gpu_data.isLumped) {
    kernelMultiply_DAMP<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_diagM(), gpu_data.get_Klocals(), p_ptr,
                                                              gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2, gpu_data.SLAU_matrix_coefs.val3,
                                                              gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, true);
  } else {
    kernelMultiply_DAMP<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Mlocals(), gpu_data.get_Klocals(), p_ptr,
                                                              gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2, gpu_data.SLAU_matrix_coefs.val3,
                                                              gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, false);
  }
}

void gpuMultiplyAlocalByVec_DAMP(gpuDataKeeper &gpu_data, int elementsCount) {
  CheckRunTime(__func__)
  float* u_ptr = gpu_data.get_u();//thrust::raw_pointer_cast( u.data() );
  float* p_ptr = gpu_data.get_p();//thrust::raw_pointer_cast( p.data() );
  kernelMultiply<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Klocals(), p_ptr);
}

void gpuCountNAdjElem(gpuDataKeeper &gpuD, int grid_size) {
  CheckRunTime(__func__)
  thrust::device_vector<float> mask_sorted(grid_size);

  //mask_sorted = d_mask;
  thrust::copy(thrust::device_pointer_cast(gpuD.get_mask()),
               thrust::device_pointer_cast(gpuD.get_mask() + grid_size),
               mask_sorted.begin());

  thrust::sort(mask_sorted.begin(), mask_sorted.end());

  thrust::reduce_by_key(mask_sorted.begin(), mask_sorted.end(), thrust::make_constant_iterator(1.f),
                        thrust::make_discard_iterator(), thrust::device_pointer_cast(gpuD.get_n_adjelem()));
}

__global__ void kernelGenerateMask(float *elements, float *mask, int DIM, int elementsCount) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < elementsCount) {
    if (DIM == 2) {
      mask[6 * (DIM - 1) * i + 0] = elements[3 * i + 0] * DIM + 0;
      mask[6 * (DIM - 1) * i + 1] = elements[3 * i + 0] * DIM + 1;
      mask[6 * (DIM - 1) * i + 2] = elements[3 * i + 1] * DIM + 0;
      mask[6 * (DIM - 1) * i + 3] = elements[3 * i + 1] * DIM + 1;
      mask[6 * (DIM - 1) * i + 4] = elements[3 * i + 2] * DIM + 0;
      mask[6 * (DIM - 1) * i + 5] = elements[3 * i + 2] * DIM + 1;
    } else if (DIM == 3) {
      mask[6 * (DIM - 1) * i + 0] = elements[4 * i + 0] * DIM + 0;
      mask[6 * (DIM - 1) * i + 1] = elements[4 * i + 0] * DIM + 1;
      mask[6 * (DIM - 1) * i + 2] = elements[4 * i + 0] * DIM + 2;
      mask[6 * (DIM - 1) * i + 3] = elements[4 * i + 1] * DIM + 0;
      mask[6 * (DIM - 1) * i + 4] = elements[4 * i + 1] * DIM + 1;
      mask[6 * (DIM - 1) * i + 5] = elements[4 * i + 1] * DIM + 2;
      mask[6 * (DIM - 1) * i + 6] = elements[4 * i + 2] * DIM + 0;
      mask[6 * (DIM - 1) * i + 7] = elements[4 * i + 2] * DIM + 1;
      mask[6 * (DIM - 1) * i + 8] = elements[4 * i + 2] * DIM + 2;
      mask[6 * (DIM - 1) * i + 9] = elements[4 * i + 3] * DIM + 0;
      mask[6 * (DIM - 1) * i + 10] = elements[4 * i + 3] * DIM + 1;
      mask[6 * (DIM - 1) * i + 11] = elements[4 * i + 3] * DIM + 2;
    }
  }
}

void gpuGenerateMask(gpuDataKeeper &gpuD, int DIM, int elementsCount) {
  CheckRunTime(__func__)
  float* elements = gpuD.get_Elements(); //thrust::raw_pointer_cast( CudaElements.data() );
  float* mask = gpuD.get_mask();//thrust::raw_pointer_cast( mask_vec.data() );
  kernelGenerateMask<<<(elementsCount + 255) / 256, 256>>>(elements, mask, DIM, elementsCount);
}

// needed for implicit scheme (isLumped = false)
__global__ void kernelCalculateDiag_DAMP(int elementsCount, float *diag, float *Mlocals, float *Klocals,
                                         float cM, float cK, float cC, float damping_alpha, float damping_beta,
                                         bool isLumped) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    for (int j = 0; j < 6; ++j) {
      if (isLumped) {
        diag[j + index * 6] =
            cM*Mlocals[j + index * 6] + cK*Klocals[j + j * 6 + index * 36] +
            cC*(damping_alpha*Mlocals[j + index * 6] + damping_beta*Klocals[j + j * 6 + index * 36]);
      } else {
        diag[j + index * 6] =
            cM*Mlocals[j + j * 6 + index * 36] + cK*Klocals[j + j * 6 + index * 36] +
            cC*(damping_alpha*Mlocals[j + j * 6 + index * 36] + damping_beta*Klocals[j + j * 6 + index * 36]);
      }
    }
  }
}

void gpuCalculateDiag_DAMP(gpuDataKeeper_DYN_DAMP &gpu_data, int elementsCount) {
  if (gpu_data.isLumped) {
    kernelCalculateDiag_DAMP<<<(elementsCount + 255) / 256, 256>>>(elementsCount,
                                                              gpu_data.get_diag(), gpu_data.get_diagM(), gpu_data.get_Klocals(),
                                                              gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2, gpu_data.SLAU_matrix_coefs.val3,
                                                              gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, true);
  } else {
    kernelCalculateDiag_DAMP<<<(elementsCount + 255) / 256, 256>>>(elementsCount,
                                                              gpu_data.get_diag(), gpu_data.get_Mlocals(), gpu_data.get_Klocals(),
                                                              gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2, gpu_data.SLAU_matrix_coefs.val3,
                                                              gpu_data.damping_coefs.val1, gpu_data.damping_coefs.val2, false);
  }
}

// needed for implicit scheme (isLumped = false)
__global__ void kernelCalculateDiag(int elementsCount, float *diag, float *Mlocals, float *Klocals,
                                    float cM, float cK) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    for (int j = 0; j < 6; ++j) {
      diag[j + index * 6] =
          cM*Mlocals[j + j * 6 + index * 36] + cK*Klocals[j + j * 6 + index * 36];
    }
  }
}


void gpuCalculateDiag(gpuDataKeeper_DYN &gpu_data, int elementsCount) {
  kernelCalculateDiag<<<(elementsCount + 255) / 256, 256>>>(elementsCount,
                                                            gpu_data.get_diag(), gpu_data.get_Mlocals(), gpu_data.get_Klocals(),
                                                            gpu_data.SLAU_matrix_coefs.val1, gpu_data.SLAU_matrix_coefs.val2);
}

void gpuDataKeeper::copyElementsFromHost(thrust::host_vector<float> &v) {
  CheckRunTime(__func__)
  this->gpuElements = v;
//  cudaMemcpy(this->get_Elements(), thrust::raw_pointer_cast(v.data()), v.size() * sizeof(float), cudaMemcpyHostToDevice);
}

void gpuDataKeeper::copyFlocalFromHost(thrust::host_vector<float> &v) {
  CheckRunTime(__func__)
  this->r = v;
}

void gpuDataKeeper::copyLoadsFromHost(thrust::host_vector<float> &v) {
  CheckRunTime(__func__)
  this->loads = v;
//  cudaMemcpy(this->get_loads(), thrust::raw_pointer_cast(v.data()), v.size() * sizeof(float), cudaMemcpyHostToDevice);

}

void gpuCopyDeviceToDevice(float *data, float *dest, int size) {
  thrust::copy(thrust::cuda::par, thrust::device_pointer_cast(data),
               thrust::device_pointer_cast(data + size), thrust::device_pointer_cast(dest));
}

void gpuCopyDeviceToHost(float *data, float *dest, int size) {
  thrust::copy(thrust::device_pointer_cast(data), thrust::device_pointer_cast(data + size), dest);
}

void gpuCopyHostToDevice(float *data, float *dest, int size) {
  thrust::copy(data, data + size, thrust::device_pointer_cast(dest));
}

void copyElementsAndFlocals(FEMdataKeeper &FEMdata, gpuDataKeeper &gpuD) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  thrust::host_vector<float> HostElements((DIM + 1) * FEMdata.elementsCount);
  thrust::host_vector<float> HostFlocal(6 * (DIM - 1) * FEMdata.elementsCount);

  int k = 0;
  for (int i = 0; i < (DIM + 1) * FEMdata.elementsCount - DIM; i += (DIM + 1)) {
    HostElements[i + 0] = FEMdata.elements[k].nodesIds[0];
    HostElements[i + 1] = FEMdata.elements[k].nodesIds[1];
    HostElements[i + 2] = FEMdata.elements[k].nodesIds[2];
    if (DIM == 3)
      HostElements[i + 3] = FEMdata.elements[k].nodesIds[3];

//    FEMdata.elements[k].Flocal.Show();

    HostFlocal[6 * (DIM - 1) * k + 0] = FEMdata.elements[k].Flocal[0];
    HostFlocal[6 * (DIM - 1) * k + 1] = FEMdata.elements[k].Flocal[1];
    HostFlocal[6 * (DIM - 1) * k + 2] = FEMdata.elements[k].Flocal[2];
    HostFlocal[6 * (DIM - 1) * k + 3] = FEMdata.elements[k].Flocal[3];
    HostFlocal[6 * (DIM - 1) * k + 4] = FEMdata.elements[k].Flocal[4];
    HostFlocal[6 * (DIM - 1) * k + 5] = FEMdata.elements[k].Flocal[5];

    if (DIM == 3) {
      HostFlocal[6 * (DIM - 1) * k + 6] = FEMdata.elements[k].Flocal[6];
      HostFlocal[6 * (DIM - 1) * k + 7] = FEMdata.elements[k].Flocal[7];
      HostFlocal[6 * (DIM - 1) * k + 8] = FEMdata.elements[k].Flocal[8];
      HostFlocal[6 * (DIM - 1) * k + 9] = FEMdata.elements[k].Flocal[9];
      HostFlocal[6 * (DIM - 1) * k + 10] = FEMdata.elements[k].Flocal[10];
      HostFlocal[6 * (DIM - 1) * k + 11] = FEMdata.elements[k].Flocal[11];
    }

    k++;
  }

  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostElements.data()), gpuD.get_Elements(), (DIM + 1) * FEMdata.elementsCount);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostFlocal.data()), gpuD.get_Flocals(), 6 * (DIM - 1) * FEMdata.elementsCount);
}

void copyFlocals(FEMdataKeeper &FEMdata, gpuDataKeeper &gpuD) {
  CheckRunTime(__func__)
  int DIM = FEMdata.DIM;
  thrust::host_vector<float> HostFlocal(3 * DIM * FEMdata.elementsCount);

  for (int k = 0; k < FEMdata.elementsCount; ++k) {
    HostFlocal[3 * DIM * k + 0] = FEMdata.elements[k].Flocal[0];
    HostFlocal[3 * DIM * k + 1] = FEMdata.elements[k].Flocal[1];
    HostFlocal[3 * DIM * k + 2] = FEMdata.elements[k].Flocal[2];
    HostFlocal[3 * DIM * k + 3] = FEMdata.elements[k].Flocal[3];
    HostFlocal[3 * DIM * k + 4] = FEMdata.elements[k].Flocal[4];
    HostFlocal[3 * DIM * k + 5] = FEMdata.elements[k].Flocal[5];
    //std::cout << HostFlocal[3 * DIM * k + 0] << " ";
  }
  //std::cout << "\n";

//  gpuD.copyFlocalFromHost(HostFlocal);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostFlocal.data()), gpuD.get_Flocals(), 6 * FEMdata.elementsCount);
}

void copyLoads(gpuDataKeeper &gpuD, std::unordered_map <int, CPU_Matrix> &loadVectors, int DIM, int elementsCount) {
  CheckRunTime(__func__)
  thrust::host_vector<float> hLoad(6 * elementsCount, 0.0f);

  for (auto& it: loadVectors) {
      hLoad[3 * DIM * it.first + 0] = loadVectors[it.first][0];
      hLoad[3 * DIM * it.first + 1] = loadVectors[it.first][1];
      hLoad[3 * DIM * it.first + 2] = loadVectors[it.first][2];
      hLoad[3 * DIM * it.first + 3] = loadVectors[it.first][3];
      hLoad[3 * DIM * it.first + 4] = loadVectors[it.first][4];
      hLoad[3 * DIM * it.first + 5] = loadVectors[it.first][5];
  }

  gpuD.copyLoadsFromHost(hLoad);
  //gpuCopyHostToDevice(thrust::raw_pointer_cast(hLoad.data()), gpuD.get_loads(), 6 * elementsCount);
}

struct AddWeighted_functor
{
    const float f1, f2;

    AddWeighted_functor(float _f1, float _f2) : f1(_f1), f2(_f2) {}

    __host__ __device__
        float operator()(const float& x, const float& y) const {
            return f1 * x + f2 * y;
        }
};

void thrustAddWeighted(thrust::device_vector<float> &v1, thrust::device_vector<float> &v2, float f1, float f2) {
    // https://stackoverflow.com/questions/22740694/multiply-device-vector-by-constant
    CheckRunTime(__func__)

    thrust::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), AddWeighted_functor(f1, f2));
}

// https://thrust.github.io/doc/group__transformed__reductions_ga0d4232a9685675f488c3cc847111e48d.html
template<typename T>
struct absolute_value : public unary_function<T,T>
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x < T(0) ? -x : x;
  }
};

float gpuCNorm(float *v, int size) {
  return thrust::transform_reduce(thrust::device_pointer_cast(v), thrust::device_pointer_cast(v + size),
                                        absolute_value<float>(),
                                        0.f,
                                        thrust::maximum<float>());
}

float thrustCNorm(thrust::device_vector<float> &v) {
    return thrust::transform_reduce(v.begin(), v.end(),
                                          absolute_value<float>(),
                                          0.f,
                                          thrust::maximum<float>());
}

void thrustGenerateMask(FEMdataKeeper FEMdata, thrust::device_vector<float> &mask) {
  CheckRunTime(__func__)
  // Make parallel!
  //thrust::device_ptr ptr_mask(mask.size());
  int DIM = FEMdata.DIM;
  thrust::host_vector<float> ptr_mask(mask.size());
  for (int eIdx = 0; eIdx < FEMdata.elementsCount; ++eIdx) {
      ptr_mask[3 * DIM * eIdx + 0] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[0] * DIM + 0);
      ptr_mask[3 * DIM * eIdx + 1] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[0] * DIM + 1);
      ptr_mask[3 * DIM * eIdx + 2] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[1] * DIM + 0);
      ptr_mask[3 * DIM * eIdx + 3] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[1] * DIM + 1);
      ptr_mask[3 * DIM * eIdx + 4] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[2] * DIM + 0);
      ptr_mask[3 * DIM * eIdx + 5] = static_cast<float>(FEMdata.elements[eIdx].nodesIds[2] * DIM + 1);
  }

  mask = ptr_mask;
}

void gpuSolveDiag(float *diag, float *r, float *res,
                  float *mask, float *n_adjelem,
                  int grid_size, int n_gl_dofs,
                  bool doAssemblyRes) {
  CheckRunTime(__func__)

  thrust::device_vector<float> diag_assemblied(n_gl_dofs), r_assemblied(n_gl_dofs);

  gpuReductionWithMask2(diag, mask, grid_size, thrust::raw_pointer_cast(diag_assemblied.data()));
  gpuReductionWithMask2(r, mask, grid_size, thrust::raw_pointer_cast(r_assemblied.data()));

  if (doAssemblyRes) {
    gpuDivide_res(thrust::raw_pointer_cast(r_assemblied.data()),
                  thrust::raw_pointer_cast(diag_assemblied.data()),
                  n_gl_dofs,
                  res);       // res = r_assemblied./diag_assemblied
  } else {    // 2N -> 6E
    thrust::device_vector<float> temp(n_gl_dofs);
    gpuDivide_res(thrust::raw_pointer_cast(r_assemblied.data()),
                  thrust::raw_pointer_cast(diag_assemblied.data()),
                  n_gl_dofs,
                  thrust::raw_pointer_cast(temp.data()));       // temp = r_assemblied./diag_assemblied
    //gpuTransform_2N_to_6E(thrust::raw_pointer_cast(temp.data()), n_gl_dofs, mask, res, grid_size);
    gpuTransformFrom_2N_to_6E(mask, thrust::raw_pointer_cast(temp.data()), res, grid_size);
  }
}

void thrustSolveDiag(thrust::device_vector<float> &diag,
                     thrust::device_vector<float> &r,
                     thrust::device_vector<float> &res,
                     thrust::device_vector<float> &mask,
                     thrust::device_vector<float> &n_adjelem,
                     bool doAssemblyRes) {
    CheckRunTime(__func__)
    int grid_size  = diag.size();
    int n_gl_dofs = n_adjelem.size();

    thrust::device_vector<float> diag_assemblied(n_gl_dofs), r_assemblied(n_gl_dofs);

    thrustReductionWithMask(diag, mask, diag_assemblied);
    thrustReductionWithMask(r, mask, r_assemblied);

    if (doAssemblyRes) {
        thrust::transform(thrust::cuda::par,
                          r_assemblied.begin(),
                          r_assemblied.end(),
                          diag_assemblied.begin(),
                          res.begin(),
                          thrust::divides<float>());         // res = r_assemblied./diag_assemblied
    } else {    // 2N -> 6E
        thrust::device_vector<float> temp(n_gl_dofs);
        thrust::transform(thrust::cuda::par,
                          r_assemblied.begin(),
                          r_assemblied.end(),
                          diag_assemblied.begin(),
                          temp.begin(),
                          thrust::divides<float>());         // temp = r_assemblied./diag_assemblied

        thrustTransform_2N_to_6E(temp, mask, res);
    }
}

void TEST_THRUST_thrustReductionWithMask() {
    const int N = 18, M = 10; // 6E
    float v[N]      = {1.f, 3.f, 3.f, 3.f, 2.f, 2.f,
                       1.f, 7.f, -5.f, -1.f, 8.f, 2.f,
                      -3.f, 1.f, 1.f, 3.f, -9.f, 1.f};
    float mask[N]   = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f,
                       5.f, 6.f, 3.f, 4.f, 7.f, 8.f,
                       7.f, 8.f, 9.f, 10.f, 5.f, 6.f};

    thrust::device_vector<float> d_v(v, v + N);
    thrust::device_vector<float> d_mask(mask, mask + N);
    thrust::device_vector<float> res(M);

    thrustReductionWithMask(d_v, d_mask, res);

    for (int i = 0; i < M; ++i) {
        float tmp = res[i];
        printf("%f ", tmp);
    }
    printf("\n");
}

// same exaple as in TEST_THRUST_thrustReductionWithMask()
void TEST_THRUST_thrustTransform_2N_to_6E() {
    const int N = 18, M = 10; // 6E
    float mask[N]   = {0.f, 1.f, 2.f, 3.f, 4.f, 5.f,
                       4.f, 5.f, 2.f, 3.f, 6.f, 7.f,
                       6.f, 7.f, 8.f, 9.f, 4.f, 5.f};
    float v[M]      = {1.f, 3.f, -2.f, 2.f, -6.f, 10.f, 5.f, 3.f, 1.f, 3.f};

    thrust::device_vector<float> d_v(v, v + M);
    thrust::device_vector<float> d_mask(mask, mask + N);
    thrust::device_vector<float> res(N);

    thrustTransform_2N_to_6E(d_v, d_mask, res);

    for (int i = 0; i < N; ++i) {
        float tmp = res[i];
        printf("%f ", tmp);
    }
    printf("\n");
}

void TEST_THRUST_CNorm() {
    const int N = 9;
    float A[N] = {1.f, 3.f, -4.f, -3.f, -2.f, 2.f, 1.f, 8.f, -9.f};

    thrust::device_vector<float> d_A(A, A + N);

    printf("%f\n", thrustCNorm(d_A));
}
