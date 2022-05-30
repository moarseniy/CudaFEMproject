#include "init.h"
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

#include <iostream>
#include <cuda.h>
#include <cusolverSp.h>

#include <math.h>
#include <vector>
#include "init.h"

#include <thrust/reduce.h>

#include "Linal2.h"
// ------------------ for gpuCalculateFEM_dyn_explicit ------------------------
#include "datakeeper.h"
#include "femfunc.h"
#include "femstruct.h"
// ----------------------------------------------------------------------------

using namespace std;

#define CUDA_CHECK_ERROR(err) \
  if (err != cudaSuccess) { \
  printf("Cuda error: %s\n", cudaGetErrorString(err)); \
  printf("Error in file: %s, line: %i\n", __FILE__, __LINE__); \
} \


gpuDataKeeper::gpuDataKeeper(int elementsCount, int nodesCount, bool doAssemblyRes) : diag(3 * DIM * elementsCount), r(3 * DIM * elementsCount),
                                                      m(3 * DIM * elementsCount), z(3 * DIM * elementsCount), s(3 * DIM * elementsCount),
                                                      p(3 * DIM * elementsCount), u(3 * DIM * elementsCount),
                                                      x(3 * DIM * elementsCount, 0.0f), mask(3 * DIM * elementsCount),
                                                      n_adjelem(DIM * nodesCount), gpuB(3 * 6 * elementsCount, 0.0f), gpuElements(3 * elementsCount),
                                                      gpuKlocals(6 * 6 * elementsCount, 0.0f), gpuFlocals(3 * DIM * elementsCount, 0.0f), tmp(3 * DIM * elementsCount),
                                                      loads(3 * DIM * elementsCount, 0.0f)
{
  CheckRunTime(__func__)
  if (doAssemblyRes)
    temp_res.resize(DIM * nodesCount);
}

gpuDataKeeper::~gpuDataKeeper() {}

//WeightedAddCoef::WeightedAddCoef(float v1, float v2): val1(v1), val2(v2) {}
//WeightedAddCoef::WeightedAddCoef(float v1, float v2, float v3): val1(v1), val2(v2), val3(v3) {}

gpuDataKeeper_DYN::gpuDataKeeper_DYN(int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped) :
  gpuDataKeeper( elementsCount, nodesCount, doAssemblyRes ),
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

gpuDataKeeper_DYN_DAMP::gpuDataKeeper_DYN_DAMP(int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped, float damping_alpha, float damping_beta) :
  gpuDataKeeper_DYN( elementsCount, nodesCount, doAssemblyRes, isLumped ) {
  this->damping_coefs.val1 = damping_alpha;
  this->damping_coefs.val2 = damping_beta;
}

gpuDataKeeper_DYN_DAMP::~gpuDataKeeper_DYN_DAMP() {}

struct divideElementwise {
  float *gpu_data;
  float value;
  divideElementwise(float *s, float v){gpu_data=s; value=v;}
  __host__ __device__
  void operator()(int i)
  {
    gpu_data[i] /= value;
  }
};

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
float det3x3(float *c) {
  return c[0]*c[4]*c[8] +
      c[1]*c[6]*c[5] +
      c[2]*c[3]*c[7] -
      c[6]*c[4]*c[2] -
      c[0]*c[5]*c[7] -
      c[1]*c[3]*c[8];
}

bool Difference_MyArray_Thrust(MyArray a, thrust::device_vector<float> b) {
    MyArray _b(a.get_size());
    for (int i = 0; i < _b.get_size(); ++i) {
        float tmp = b[i];
        _b[i] = tmp;
    }

    return a.equalsToArray(_b, 1e-15);
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
}

#define id(m, i, j) m[j + 3 * i]

__device__
void inverse3x3(float *ic, float *c, float det) {
  float invdet = 1.0f / det;
  ic[0 + 3 * 0] = (c[1 + 3 * 1] * c[2 + 3 * 2] - c[1 + 3 * 2] * c[2 + 3 * 1]) * invdet;
  ic[1 + 3 * 0] = (c[2 + 3 * 0] * c[1 + 3 * 2] - c[1 + 3 * 0] * c[2 + 3 * 2]) * invdet;
  ic[2 + 3 * 0] = (c[1 + 3 * 0] * c[2 + 3 * 1] - c[2 + 3 * 0] * c[1 + 3 * 1]) * invdet;

  ic[0 + 3 * 1] = (c[2 + 3 * 1] * c[0 + 3 * 2] - c[0 + 3 * 1] * c[2 + 3 * 2]) * invdet;
  ic[1 + 3 * 1] = (c[0 + 3 * 0] * c[2 + 3 * 2] - c[2 + 3 * 0] * c[0 + 3 * 2]) * invdet;
  ic[2 + 3 * 1] = (c[0 + 3 * 1] * c[2 + 3 * 0] - c[0 + 3 * 0] * c[2 + 3 * 1]) * invdet;

  ic[0 + 3 * 2] = (c[0 + 3 * 1] * c[1 + 3 * 2] - c[0 + 3 * 2] * c[1 + 3 * 1]) * invdet;
  ic[1 + 3 * 2] = (c[0 + 3 * 2] * c[1 + 3 * 0] - c[0 + 3 * 0] * c[1 + 3 * 2]) * invdet;
  ic[2 + 3 * 2] = (c[0 + 3 * 0] * c[1 + 3 * 1] - c[0 + 3 * 1] * c[1 + 3 * 0]) * invdet;
}

__device__
float det3(float a0, float a1, float a2, float a3,
           float a4, float a5, float a6, float a7, float a8) {
  return a0*a4*a8 +
      a1*a6*a5 +
      a2*a3*a7 -
      a6*a4*a2 -
      a0*a5*a7 -
      a1*a3*a8;
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
    float v1 = det3(c[5], c[6], c[7], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v2 = det3(c[1], c[2], c[3], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v3 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[13], c[14], c[15]);
    float v4 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[9], c[10], c[11]);
    return v1 - v2 + v3 - v4;
  }
}

__device__
void Get_matrix(float *a, int n, float *c, int indRow, int indCol) {
  //float *a = (float*)malloc(3 * 3 * sizeof (float));
  int ki = 0;
  for (int i = 0; i < n; i++) {
    if (i != indRow) {
      for (int j = 0, kj = 0; j < n; j++) {
        if (j != indCol) {
          a[kj + ki * 3] = c[j + i * n];
          kj++;
        }
      }
      ki++;
    }
  }

  //return a;
}

__device__
void inverse(float *ic, float *b, int size) {
  //float *ic = (float*)malloc(4 * 4 * sizeof (float));
  //printf("%f", b[0]);
  float determinant = det3x3(b);

  if (determinant) {
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        float *temp =(float*)malloc(3 * 3 * sizeof (float));
        //__shared__ float temp[3 * 3];
        Get_matrix(temp, size, b, i, j);
        ic[j + i * size] = ((i + j + 2) % 2 == 0 ? 1.0 : -1.0) * det(temp, 2) / determinant;
        free(temp);
      }
    }
  }

  float swap;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (i > j) {
        swap = ic[j + i * size];
        ic[j + i * size] = ic[i + j * size];
        ic[i + j * size] = swap;
      }
    }
  }

  //return ic;
}

//void TEST_THRUST() {
//  const int N = 7;
//  int A[N] = {1, 3, 3, 3, 2, 2, 1}; // input keys
//  int B[N] = {9, 8, 7, 6, 5, 4, 3}; // input values
//  int C[N];                         // output keys
//  int D[N];                         // output values
//  thrust::pair<int*,int*> new_end;
//  thrust::equal_to<int> binary_pred;
//  thrust::plus<int> binary_op;
//  new_end = thrust::reduce_by_key(A, A + N, B, C, D, binary_pred, binary_op);
//  for (int i = 0; i < 4; ++i) {
//    std::cout << C[i] << " ";
//  }
//  std::cout << std::endl;
//  for (int i = 0; i < 4; ++i) {
//    std::cout << D[i] << " ";
//  }
//}

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
  //    temp_mask = mask;
  //    temp_v = v;

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
        host_res[i] = host_v[static_cast<int>(host_mask[i])];
    }

    gpuCopyHostToDevice(thrust::raw_pointer_cast(host_res.data()), d_res, grid_size);
//    gpuCopy(thrust::raw_pointer_cast(host_res.data()), d_res, grid_size);

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
  //    thrust::transform(a.begin(), a.end(), b.begin(), c.begin(),
  //                      thrust::placeholders::_1 / thrust::placeholders::_2);
  thrust::transform(thrust::cuda::par,
                    thrust::device_pointer_cast(a),
                    thrust::device_pointer_cast(a + size),
                    thrust::device_pointer_cast(b),
                    thrust::device_pointer_cast(c),
                    thrust::divides<float>());
}

__global__ void kernelCalculateKlocal2(int elementsCount, float *elements, float *diag,
                                       float *nodesX, float *nodesY, int nodesCount,
                                       float *D, float *K, float *B,
                                       float *constraints, int constraintsCount) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < elementsCount) {
    float *C = (float*)malloc(3 * 3 * sizeof(float));
    float *IC = (float*)malloc(3 * 3 * sizeof(float));

    C[0 + 0 * 3] = C[0 + 1 * 3] = C[0 + 2 * 3] = 1.0f;
    C[1 + 0 * 3] = nodesX[int(elements[3 * index + 0])];
    C[1 + 1 * 3] = nodesX[int(elements[3 * index + 1])];
    C[1 + 2 * 3] = nodesX[int(elements[3 * index + 2])];
    C[2 + 0 * 3] = nodesY[int(elements[3 * index + 0])];
    C[2 + 1 * 3] = nodesY[int(elements[3 * index + 1])];
    C[2 + 2 * 3] = nodesY[int(elements[3 * index + 2])];

    float determinant = det3x3(C);

    inverse3x3(IC, C, determinant);
    //inverse(IC, C, 3);

    //float *B = (float*)malloc(3 * 6 * sizeof(float));
    for (int i = 0; i < 3; i++) {
      B[(2 * i + 0) + 0 * 6 + index * 3 * 6] = IC[i + 3 * 1];
      B[(2 * i + 1) + 0 * 6 + index * 3 * 6] = 0.0f;
      B[(2 * i + 0) + 1 * 6 + index * 3 * 6] = 0.0f;
      B[(2 * i + 1) + 1 * 6 + index * 3 * 6] = IC[i + 3 * 2];
      B[(2 * i + 0) + 2 * 6 + index * 3 * 6] = IC[i + 3 * 2];
      B[(2 * i + 1) + 2 * 6 + index * 3 * 6] = IC[i + 3 * 1];
    }



    float *B_T = (float*)malloc(6 * 3 * sizeof(float));
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++) {
        B_T[j + i * 3] = B[i + j * 6 + index * 3 * 6];
      }
    }

    float *temp = (float*)malloc(6 * 3 * sizeof(float));
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++) {
        temp[j + i * 3] = 0.0f;
        for (int k = 0; k < 3; k++) {
          temp[j + i * 3] += B_T[k + i * 3] * D[j + k * 3];
        }
      }
    }

    //float *K = (float*)malloc(6 * 6 * sizeof(float));
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        K[j + i * 6 + 36 * index] = 0.0f;
        for (int k = 0; k < 3; k++) {
          K[j + i * 6 + 36 * index] += temp[k + i * 3] * B[j + k * 6 + index * 3 * 6];
        }
      }
    }

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

    //      printf("GPU(p)%f \n", p[5]);
    //      //printf("CUDA\n");
    //      for (int j = 0; j < 6; j++) {
    //        u[j + index * 6] = 0.0f;
    //        for (int k = 0; k < 6; k++) {
    //            u[j + index * 6] += K[k + j * 6] * p[k + index * 6];
    //            if (j + index * 6 == 0)
    //                printf("%f += %f * %f\n", u[0], K[k + j * 6], p[k + index * 6]);
    //        }
    //      }
  }

}

void gpuCalculateKlocal2(gpuDataKeeper &gpu_data, int elementsCount,
                         float *h_nodesX, float *h_nodesY, int nodesCount,
                         float *h_D, float *h_constraints, int constraintsCount) {
  CheckRunTime(__func__)

  float* elements = gpu_data.get_Elements();//thrust::raw_pointer_cast( CudaElements.data() );
  float* diag = gpu_data.get_diag();//thrust::raw_pointer_cast( diag_vec.data() );

  float *d_elements, *d_constraints;
  float *d_nodesX, *d_nodesY, *d_p, *d_u, *d_D;

  int grid_size = 6 * elementsCount;

  cudaMalloc((void**)&d_D, 3 * 3 * sizeof(float));
  cudaMemcpy(d_D, h_D, 3 * 3 * sizeof(float), cudaMemcpyHostToDevice);

  //cudaMalloc((void**)&d_elements, 3 * elementsCount * sizeof(float));
  //cudaMemcpy(d_elements, h_elements, 3 * elementsCount * sizeof(float), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_constraints, constraintsCount * sizeof(float));
  cudaMemcpy(d_constraints, h_constraints, constraintsCount * sizeof(float), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&d_nodesX, nodesCount * sizeof(float));
  cudaMalloc((void**)&d_nodesY, nodesCount * sizeof(float));
  cudaMemcpy(d_nodesX, h_nodesX, nodesCount * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nodesY, h_nodesY, nodesCount * sizeof(float), cudaMemcpyHostToDevice);

  //    cudaMalloc((void**)&d_p, grid_size * sizeof(float));
  //    cudaMalloc((void**)&d_u, grid_size * sizeof(float));
  //    cudaMemcpy(d_p, h_p, grid_size * sizeof(float), cudaMemcpyHostToDevice);
  //    cudaMemcpy(d_u, h_u, grid_size * sizeof(float), cudaMemcpyHostToDevice);

  kernelCalculateKlocal2<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements, diag,
                                                               d_nodesX, d_nodesY, nodesCount,
                                                               d_D, gpu_data.get_Klocals(), gpu_data.get_B(),
                                                               d_constraints, constraintsCount);
  //    cudaFree(d_p);
  //    cudaFree(d_u);
  cudaFree(d_nodesX);
  cudaFree(d_nodesY);
  //    cudaFree(d_elements);
  cudaFree(d_constraints);
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

void gpuCalculateMlocal(gpuDataKeeper_DYN &gpu_data, int elementsCount,
                         float *h_nodesX, float *h_nodesY, int nodesCount, float rho) {
  CheckRunTime(__func__)
  float* elements = gpu_data.get_Elements();
  float *d_nodesX, *d_nodesY;

  cudaMalloc((void**)&d_nodesX, nodesCount * sizeof(float));
  cudaMalloc((void**)&d_nodesY, nodesCount * sizeof(float));
  cudaMemcpy(d_nodesX, h_nodesX, nodesCount * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nodesY, h_nodesY, nodesCount * sizeof(float), cudaMemcpyHostToDevice);

  if (gpu_data.isLumped) {
    kernelCalculateMlocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements,
                                                                d_nodesX, d_nodesY, gpu_data.get_diagM(), rho, true);
//    kernelCalculateMlocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements,
//                                                                d_nodesX, d_nodesY, gpu_data.get_Mlocals(), rho, true);
  } else {
    kernelCalculateMlocal<<<(elementsCount + 255) / 256, 256>>>(elementsCount, elements,
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
          if (k == j) factor += damping_alpha*Mlocals[k + index * 6];
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


void gpuMultiplyKlocalByVec(gpuDataKeeper &gpu_data, int elementsCount) {
  CheckRunTime(__func__)
  float* u_ptr = gpu_data.get_u();//thrust::raw_pointer_cast( u.data() );
  float* p_ptr = gpu_data.get_p();//thrust::raw_pointer_cast( p.data() );
  kernelMultiply<<<(elementsCount + 255) / 256, 256>>>(elementsCount, u_ptr, gpu_data.get_Klocals(), p_ptr);
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
  thrust::copy(thrust::device, thrust::device_pointer_cast(gpuD.get_mask()),
               thrust::device_pointer_cast(gpuD.get_mask() + grid_size),
               mask_sorted.begin());

  thrust::sort(mask_sorted.begin(), mask_sorted.end());

  thrust::reduce_by_key(mask_sorted.begin(), mask_sorted.end(), thrust::make_constant_iterator(1.f),
                        thrust::make_discard_iterator(), thrust::device_pointer_cast(gpuD.get_n_adjelem()));
}

__global__ void kernelGenerateMask(float *elements, float *mask, int elementsCount) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < elementsCount) {
    mask[3 * DIM * i + 0] = elements[3 * i + 0] * DIM + 0;
    mask[3 * DIM * i + 1] = elements[3 * i + 0] * DIM + 1;
    mask[3 * DIM * i + 2] = elements[3 * i + 1] * DIM + 0;
    mask[3 * DIM * i + 3] = elements[3 * i + 1] * DIM + 1;
    mask[3 * DIM * i + 4] = elements[3 * i + 2] * DIM + 0;
    mask[3 * DIM * i + 5] = elements[3 * i + 2] * DIM + 1;
  }
}

void gpuGenerateMask(gpuDataKeeper &gpuD, int elementsCount) {
  CheckRunTime(__func__)
  float* elements = gpuD.get_Elements(); //thrust::raw_pointer_cast( CudaElements.data() );
  float* mask = gpuD.get_mask();//thrust::raw_pointer_cast( mask_vec.data() );
  kernelGenerateMask<<<(elementsCount + 255) / 256, 256>>>(elements, mask, elementsCount);
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

void gpuDataKeeper::copyElementsFromHost(thrust::host_vector<float> v) {
  CheckRunTime(__func__)
  this->gpuElements = v;
}
void gpuDataKeeper::copyFlocalFromHost(thrust::host_vector<float> v) {
  CheckRunTime(__func__)
  this->r = v;
}

void gpuDataKeeper::setZeroVec() {
  //thrust::fill(z.begin(), z.end(), 0.0f);
  //thrust::fill(m.begin(), m.end(), 0.0f);
  //thrust::fill(p.begin(), p.end(), 0.0f);
  //thrust::fill(s.begin(), s.end(), 0.0f);
  //thrust::fill(temp_res.begin(), temp_res.end(), 0.0f);

  thrust::fill(x.begin(), x.end(), 0.0f);

  //thrust::fill(u.begin(), u.end(), 0.0f);
  //thrust::fill(r.begin(), r.end(), 0.0f);
}

void gpuDataKeeper::copyLoadsFromHost(thrust::host_vector<float> v) {
  CheckRunTime(__func__)
  this->loads = v;
}

void gpuCopyDeviceToDevice(float *data, float *dest, int size) {
  thrust::copy(thrust::device, thrust::device_pointer_cast(data),
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
  thrust::host_vector<float> HostElements(3 * FEMdata.elementsCount);
  thrust::host_vector<float> HostFlocal(6 * FEMdata.elementsCount);

  int k = 0;
  for (int i = 0; i < 3 * FEMdata.elementsCount - 2; i += 3) {
    HostElements[i + 0] = FEMdata.elements[k].nodesIds[0];
    HostElements[i + 1] = FEMdata.elements[k].nodesIds[1];
    HostElements[i + 2] = FEMdata.elements[k].nodesIds[2];

    HostFlocal[3 * DIM * k + 0] = FEMdata.elements[k].Flocal[0];
    HostFlocal[3 * DIM * k + 1] = FEMdata.elements[k].Flocal[1];
    HostFlocal[3 * DIM * k + 2] = FEMdata.elements[k].Flocal[2];
    HostFlocal[3 * DIM * k + 3] = FEMdata.elements[k].Flocal[3];
    HostFlocal[3 * DIM * k + 4] = FEMdata.elements[k].Flocal[4];
    HostFlocal[3 * DIM * k + 5] = FEMdata.elements[k].Flocal[5];

    k++;
  }
  //gpuD.copyElementsFromHost(HostElements);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostElements.data()), gpuD.get_Elements(), 3 * FEMdata.elementsCount);
  //gpuD.copyFlocalFromHost(HostFlocal);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostFlocal.data()), gpuD.get_Flocals(), 6 * FEMdata.elementsCount);
}

void copyFlocals(FEMdataKeeper &FEMdata, gpuDataKeeper &gpuD) {
  CheckRunTime(__func__)
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

  //gpuD.copyFlocalFromHost(HostFlocal);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(HostFlocal.data()), gpuD.get_Flocals(), 6 * FEMdata.elementsCount);
}

void copyLoads(gpuDataKeeper &gpuD, std::unordered_map <int, MyArray> &loadVectors, int elementsCount) {
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

//  gpuD.copyLoadsFromHost(hLoad);
  gpuCopyHostToDevice(thrust::raw_pointer_cast(hLoad.data()), gpuD.get_loads(), 6 * elementsCount);
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
    gpuTransform_2N_to_6E(thrust::raw_pointer_cast(temp.data()), n_gl_dofs, mask, res, grid_size);
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

// Rewrite it in the form that SparseMatrixCOO is stored on device, not host!
// Make parallel!
void thrustMultiplyByVector(SparseMatrixCOO &M,
                            thrust::device_vector<float> &vec,
                            thrust::device_vector<float> &res) {
    CheckRunTime(__func__)
    int *x = M.get_x(), *y = M.get_y();
    float *data = M.get_data();

    thrust::host_vector<float> host_vec(vec.size()), host_res(res.size(), 0.0f);
    host_vec = vec;

    for (int i = 0; i < M.get_size(); ++i) {
        host_res[x[i]] += data[i] * host_vec[y[i]];
    }

    res = host_res;
}

void gpuCalculateFEM_dyn_relaxation(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float dt, float beta1, float eps) {
    // Zienkiewicz, Taylor, Zhu "The Finite Element Method: Its Basis and Fundamentals" 6th edition 17.3.3 GN22 (page 608)
    CheckRunTime(__func__)

    int n_elems = FEMdata.elementsCount;
    int n_gl_dofs = FEMdata.nodesCount * DIM;
    int grid_size = 3 * DIM * n_elems;

    bool is_damping = !((damping_alpha == 0.0f) && (damping_beta == 0.0f));

    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
        Element *elem = &FEMdata.elements[eIdx];
        elem->CalculateKlocal(FEMdata.D, FEMdata.nodesX, FEMdata.nodesY);
        elem->CalculateMlocal(rho, FEMdata.nodesX, FEMdata.nodesY, true); // Lumping on
        if (is_damping) elem->CalculateClocal(damping_alpha, damping_beta);
    }

    std::unordered_map <int, std::vector<int>> nodeAdjElem;
    CalculateNodeAdjElem(FEMdata, nodeAdjElem);

    ApplyConstraints_EbE(FEMdata);
    AssignLoadElement(FEMdata, nodeAdjElem);

    thrust::device_vector<float>    x(grid_size, 0.0f), vel(grid_size, 0.0f),
                                    res(grid_size, 0.0f), b(grid_size),
                                    F(grid_size), diag(grid_size),
                                    mask(grid_size);
    thrust::device_vector<float>    n_adjelem(n_gl_dofs);

    thrust::host_vector<float>      hF(grid_size, 0.0f), hLoad(grid_size, 0.0f),
                                    hdiag(grid_size);
    thrust::device_vector<float>    dLoad(grid_size);

    thrustGenerateMask(FEMdata, mask);
    thrustCountNAdjElem(mask, n_adjelem);

    int nonzeroK = 0, nonzeroC = 0;
    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
        nonzeroK += FEMdata.elements[eIdx].Klocal.CountNonzero();
    }
    if (is_damping) {
        for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
            nonzeroC += FEMdata.elements[eIdx].Clocal.CountNonzero();
        }
    }

    SparseMatrixCOO K(nonzeroK), C(nonzeroC);
    for (int eIdx = 0; eIdx < n_elems; ++eIdx) {
        Element elem = FEMdata.elements[eIdx];

        for (int xlocal = 0; xlocal < 3 * DIM; ++xlocal) {
            for (int ylocal = 0; ylocal < 3 * DIM; ++ylocal) {
                K.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Klocal(xlocal, ylocal));
                if (is_damping)
                    C.write_value(3 * DIM * eIdx + xlocal, 3 * DIM * eIdx + ylocal, elem.Clocal(xlocal, ylocal));
            }
        }

        // Get diagonal elements of mass matrix
        float tmp_M = elem.Mlocal(0, 0);
        hdiag[3 * DIM * eIdx + 0] = tmp_M;
        hdiag[3 * DIM * eIdx + 1] = tmp_M;
        hdiag[3 * DIM * eIdx + 2] = tmp_M;
        hdiag[3 * DIM * eIdx + 3] = tmp_M;
        hdiag[3 * DIM * eIdx + 4] = tmp_M;
        hdiag[3 * DIM * eIdx + 5] = tmp_M;
    }

    diag = hdiag;

    int nt = 1;
    float cnorm_res, cnorm_vel;
    do {
        float t = nt*dt;
        std::cout << "======= Time iteration #" << nt << " =======" << std::endl;
        std::cout << "=========== Time " << t << " ===========" << std::endl << std::endl;

        thrustAddWeighted(x, vel, 1.0f, dt);
        thrustAddWeighted(x, res, 1.0f, 0.5f * dt*dt);
        thrustAddWeighted(vel, res, 1.0f, (1.0f - beta1) * dt);

        thrustMultiplyByVector(K, x, b);
        if (is_damping) {
            thrust::device_vector<float> tmp(grid_size);
            thrustMultiplyByVector(C, vel, tmp);
            thrust::transform(thrust::cuda::par,
                              tmp.begin(),
                              tmp.end(),
                              b.begin(),
                              b.begin(),
                              thrust::plus<float>());               // b = b + tmp
        }

        // update F
        for (int beIdx = 0; beIdx < FEMdata.boundaryEdgesCount; ++beIdx) {
            BoundaryEdge bedge = FEMdata.boundary[beIdx];
            int eIdx = bedge.adj_elem1;
            FEMdata.elements[eIdx].CalculateFlocal(bedge, FEMdata.nodesX, FEMdata.nodesY, t);
            hF[3 * DIM * eIdx + 0] = FEMdata.elements[eIdx].Flocal[0];
            hF[3 * DIM * eIdx + 1] = FEMdata.elements[eIdx].Flocal[1];
            hF[3 * DIM * eIdx + 2] = FEMdata.elements[eIdx].Flocal[2];
            hF[3 * DIM * eIdx + 3] = FEMdata.elements[eIdx].Flocal[3];
            hF[3 * DIM * eIdx + 4] = FEMdata.elements[eIdx].Flocal[4];
            hF[3 * DIM * eIdx + 5] = FEMdata.elements[eIdx].Flocal[5];
            // ToDO: ADD APPLY CONSTRAINTS FOR Flocal!
        }
        F = hF;

        thrust::transform(thrust::cuda::par,
                          F.begin(),
                          F.end(),
                          b.begin(),
                          b.begin(),
                          thrust::minus<float>());               // b = F - b

        // Think how not to use loadVectors! Too dificult!
        std::unordered_map <int, MyArray> loadVectors;
        loadVectors.clear();    // in order to GetMapElement2Loadvector, because if not cleared the values are added
                                // instead of assigned. See if-statement in for-loop in the function's body
        GetMapElement2Loadvector(FEMdata, loadVectors, t);
        for (auto& it: loadVectors) {
            hLoad[3 * DIM * it.first + 0] = loadVectors[it.first][0];
            hLoad[3 * DIM * it.first + 1] = loadVectors[it.first][1];
            hLoad[3 * DIM * it.first + 2] = loadVectors[it.first][2];
            hLoad[3 * DIM * it.first + 3] = loadVectors[it.first][3];
            hLoad[3 * DIM * it.first + 4] = loadVectors[it.first][4];
            hLoad[3 * DIM * it.first + 5] = loadVectors[it.first][5];
        }
        dLoad = hLoad;
        thrust::transform(thrust::cuda::par,
                          b.begin(),
                          b.end(),
                          dLoad.begin(),
                          b.begin(),
                          thrust::plus<float>());

        thrustSolveDiag(diag, b, res, mask, n_adjelem, false);

        thrustAddWeighted(vel, res, 1.0f, beta1 * dt);

        ++nt;
        cnorm_res = thrustCNorm(res);
        cnorm_vel = thrustCNorm(vel);
        std::cout << "C-norm res = " << cnorm_res << "\nC-norm vel = " << cnorm_vel << std::endl;
        std::cout << std::endl;
        if (nt >= 10000) break;
    } while (!((cnorm_res < eps) && (cnorm_vel < eps)));

    thrust::device_vector<float> displ(n_gl_dofs);
    thrustReductionWithMask(x, mask, displ);
    thrust::transform(thrust::cuda::par,
                      displ.begin(),
                      displ.end(),
                      n_adjelem.begin(),
                      displ.begin(),
                      thrust::divides<float>());         // displ = displ./n_adjelem

    thrust::copy(displ.begin(), displ.end(), FEMdata.displacements.get_data()); // is the last argument OK?
}

void TEST_THRUST_addWeighted() {
    const int N = 7;
    float A[N] = {1.f, 3.f, 3.f, 3.f, 2.f, 2.f, 1.f};
    float B[N] = {9.f, 8.f, 7.f, 6.f, 5.f, 4.f, 3.f};
    float f1 = 3.0;
    float f2 = 5.0;

    thrust::device_vector<float> d_A(A, A + N);
    thrust::device_vector<float> d_B(B, B + N);

    thrustAddWeighted(d_A, d_B, f1, f2);

    for (int i = 0; i < N; ++i) {
        float tmp = d_A[i];
        printf("%f ", tmp);
    }
    printf("\n");
}




void TEST_THRUST_GenerateMask(FEMdataKeeper &FEMdata) {
    int n_elems = FEMdata.elementsCount;
    int grid_size = 3 * DIM * n_elems;
    thrust::device_vector<float>    thrust_mask(grid_size, 0.0f);
    MyArray mask(grid_size);

    GenerateMask(FEMdata, mask);
    thrustGenerateMask(FEMdata, thrust_mask);

    for (int i = 0; i < n_elems; ++i) {
        int tmp = int(thrust_mask[i]);
        if (mask[i] != tmp) {
            printf("Error! %d != %d\n", mask[i], tmp);
            return;
        }
    }
    printf("All correct!\n");
}

void TEST_THRUST_MultiplyByVector() {
    const int N = 7;
    float A[N] = {1.f, 3.f, 3.f, 3.f, 2.f, 2.f, 1.f};

    thrust::device_vector<float> d_A(A, A + N);
    thrust::device_vector<float> d_B(N);

    SparseMatrixCOO M(10);
    M.write_value(0,0, 3.f);
    M.write_value(1,1, 3.f);
    M.write_value(2,2, 3.f);
    M.write_value(3,3, 3.f);
    M.write_value(4,4, 3.f);
    M.write_value(5,5, 3.f);
    M.write_value(6,6, 3.f);
    M.write_value(0,2, -1.f);
    M.write_value(2,4, 6.f);
    M.write_value(3,4, -3.f);

    thrustMultiplyByVector(M,d_A, d_B);

    for (int i = 0; i < N; ++i) {
        float tmp = d_B[i];
        printf("%f ", tmp);
    }
    printf("\n");
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
