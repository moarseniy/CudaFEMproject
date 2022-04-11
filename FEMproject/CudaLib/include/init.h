#ifndef CUDA_LIB_
#define CUDA_LIB_

#include <iostream>
//#include "femfunc.h"


void TEST_THRUST();
void gpuAddWeighted(float *a, float *b, float v1, float v2, int size);
void gpuReductionWithMaskAndTransform(float *v, float *mask, int size, float *res, int size_new);
void gpuReductionWithMask(float *v, float *mask, int size, float *res, int size_new);
float gpuDotProduct(float *a, float *b, int size);




#endif
