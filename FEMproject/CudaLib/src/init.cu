
#include "init.h"

using namespace std;

#define CUDA_CHECK_ERROR(err) \
  if (err != cudaSuccess) { \
  printf("Cuda error: %s\n", cudaGetErrorString(err)); \
  printf("Error in file: %s, line: %i\n", __FILE__, __LINE__); \
} \


gpuDataKeeper::gpuDataKeeper(int elementsCount, int nodesCount, bool doAssemblyRes) : diag(3 * DIM * elementsCount), r(3 * DIM * elementsCount),
                                                      m(3 * DIM * elementsCount), z(3 * DIM * elementsCount), s(3 * DIM * elementsCount),
                                                      p(3 * DIM * elementsCount), u(3 * DIM * elementsCount),
                                                      x(3 * DIM * elementsCount), mask(3 * DIM * elementsCount),
                                                      n_adjelem(DIM * nodesCount), gpuB(3 * 6 * elementsCount),
                                                      gpuKlocals(6 * 6 * elementsCount), gpuFlocals(6 * elementsCount){
  CheckRunTime(__func__)
  if (doAssemblyRes)
  temp_res.resize(DIM * nodesCount);
}

gpuDataKeeper::~gpuDataKeeper() {}

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

__device__
float det3x3(float *c) {
  return c[0]*c[4]*c[8] +
      c[1]*c[6]*c[5] +
      c[2]*c[3]*c[7] -
      c[6]*c[4]*c[2] -
      c[0]*c[5]*c[7] -
      c[1]*c[3]*c[8];
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

void gpuMultiplyKlocalByVec(gpuDataKeeper &gpu_data, int elementsCount) {
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

void gpuDataKeeper::copyElementsFromHost(thrust::host_vector<float> v) {
  CheckRunTime(__func__)
  this->gpuElements = v;
}
void gpuDataKeeper::copyFlocalFromHost(thrust::host_vector<float> v) {
  CheckRunTime(__func__)
  this->r = v;
}

void gpuCopyDeviceToDevice(float *data, float *dest, int size) {
  thrust::copy(thrust::device, thrust::device_pointer_cast(data),
               thrust::device_pointer_cast(data + size), thrust::device_pointer_cast(dest));
}

void gpuCopyDeviceToHost(float *data, float *dest, int size) {
  thrust::copy(thrust::device_pointer_cast(data), thrust::device_pointer_cast(data + size), dest);
}

void copyElements(FEMdataKeeper FEMdata, gpuDataKeeper &gpuD) {
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
  gpuD.copyElementsFromHost(HostElements);
  gpuD.copyFlocalFromHost(HostFlocal);
}

