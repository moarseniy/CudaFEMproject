#include "matrix_kernels.h"

//cublasSgemmStridedBatched(handle, CUBLAS_OP_T, CUBLAS_OP_N,
//                          2, 1, 2,
//                          &alpha,
//                          d_A, 2, 4,
//                          d_B, 2, 2,
//                          &beta,
//                          d_C, 2, 2,
//                          4);


__global__
void show(float* ptr, int size) {
  for(int i =0; i<size; i++)
  printf("%f\n", ptr[i]);
}

float thrust_dotProduct_Ker(float *src1, float *src2, size_t size) {
  return thrust::inner_product(thrust::device,
                               thrust::device_pointer_cast(src1),
                               thrust::device_pointer_cast(src1 + size),
                               thrust::device_pointer_cast(src2), 0.f);
}

void thrust_divideElementwise_Ker(float *src1, float *src2, float *tgt, size_t size) {
  thrust::transform(thrust::device,
                    thrust::device_pointer_cast(src1),
                    thrust::device_pointer_cast(src1 + size),
                    thrust::device_pointer_cast(src2),
                    thrust::device_pointer_cast(tgt),
                    thrust::divides<float>());
}

void thrust_sortByKey_Ker(float *keys, float *src, size_t size) {
  thrust::sort_by_key(thrust::device,
                      thrust::device_pointer_cast(keys),
                      thrust::device_pointer_cast(keys + size),
                      thrust::device_pointer_cast(src));
}

void thrust_reduceByKey_Ker(float *keys, float *src, float *tgt, size_t size) {
  thrust::reduce_by_key(thrust::device,
                        thrust::device_pointer_cast(keys),
                        thrust::device_pointer_cast(keys + size),
                        thrust::device_pointer_cast(src),
                        thrust::make_discard_iterator(),
                        thrust::device_pointer_cast(tgt));
}

void thrust_sort_Ker(float *src, size_t size) {
  try {
    thrust::stable_sort(thrust::device,
                        thrust::device_pointer_cast(src),
                        thrust::device_pointer_cast(src + size));
  } catch(thrust::system_error e) {
    std::cerr << "Error inside sort: " << e.what() << std::endl;
    EXIT_FAILURE;
  }
}

void thrust_setTo_Ker(float *src, size_t size, float v) {
  thrust::fill(thrust::device,
               thrust::device_pointer_cast(src),
               thrust::device_pointer_cast(src + size), v);
}

__global__
void kernelMultiply(const size_t n, const float *data, const float *vec, const size_t size, float *tgt) {
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

void multiplyByVec_Ker(const size_t numMatr, const float *data,
                       const float *vec_data, const size_t vecSize,
                       float *tgt) {
  kernelMultiply<<<(numMatr + 255) / 256, 256>>>(numMatr, data, vec_data, vecSize, tgt);
}

__global__
void kernelSetTo(const size_t size, float *data, const float value) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < size) {
    data[id] = value;
  }
}

void setTo_Ker(const size_t size, float *data, const float value) {
  kernelSetTo<<<(size + 255) / 256, 256>>>(size, data, value);
}

__global__
void kernelScale(const size_t size, float *data, const float value) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < size) {
    data[id] *= value;
  }
}

void scale_Ker(const size_t size, float *data, const float value) {
  kernelScale<<<(size + 255) / 256, 256>>>(size, data, value);
}

__global__
void kernelAddWeighted(const size_t size, float *data, const float *src,
                       const float alpha, const float beta) {
  size_t id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < size) {
    data[id] = alpha * data[id] + beta * src[id];
  }
}

void addWeighted_Ker(const size_t size, float *data, const float *src,
                     const float alpha, const float beta) {
  kernelAddWeighted<<<(size + 255) / 256, 256>>>(size, data, src, alpha, beta);
}

__host__ __device__ unsigned int
IDX(unsigned int i,unsigned  int j,unsigned int ld){
  return j*ld+i;
}

__device__ float
det_kernel(float *a_copy,unsigned int *n,cublasHandle_t *hdl){
  int *info = (int *)malloc(sizeof(int));info[0]=0;
  int batch=1;int *p = (int *)malloc(*n*sizeof(int));
  float **a = (float **)malloc(sizeof(float *));
  *a = a_copy;
  cublasStatus_t status=cublasSgetrfBatched(*hdl, *n, a, *n, p, info, batch);
  cudaDeviceSynchronize();
  unsigned int i1;
  float res=1;
  for(i1=0;i1<(*n);++i1)res*=a_copy[IDX(i1,i1,*n)];
  return res;
}
