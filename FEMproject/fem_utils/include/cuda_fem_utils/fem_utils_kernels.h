#ifndef FEM_UTILS_KERNELS_H_
#define FEM_UTILS_KERNELS_H_

#include <cstdio>
#include <cuda_runtime.h>
#include <cuda.h>

void genMask_ker(float *mask, float *elements, size_t dim, size_t elementsCount);

#endif /* FEM_UTILS_KERNELS_H_ */
