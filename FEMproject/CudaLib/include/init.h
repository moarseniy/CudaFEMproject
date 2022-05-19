#ifndef CUDA_LIB_
#define CUDA_LIB_

#include <iostream>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/inner_product.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>

#include <cuda.h>
#include <cusolverSp.h>

#include <math.h>
#include <vector>
#include "Linal2.h"
#include "Tools.h"
#include "femstruct.h"
#include "datakeeper.h"

class gpuDataKeeper{
public:
  gpuDataKeeper(int elementsCount, int nodesCount, bool doAssemblyRes);
  ~gpuDataKeeper();
  float* get_Klocals() { return thrust::raw_pointer_cast(gpuKlocals.data());}
  float* get_B() { return thrust::raw_pointer_cast(gpuB.data());}
  float* get_Flocals() { return thrust::raw_pointer_cast(gpuFlocals.data());}
  float* get_Elements() { return thrust::raw_pointer_cast(gpuElements.data());}
  float* get_r() { return thrust::raw_pointer_cast(r.data());}
  float* get_mask() { return thrust::raw_pointer_cast(mask.data());}
  float* get_n_adjelem() { return thrust::raw_pointer_cast(n_adjelem.data());}
  float* get_diag() { return thrust::raw_pointer_cast(diag.data());}
  float* get_m() { return thrust::raw_pointer_cast(m.data());}
  float* get_z() { return thrust::raw_pointer_cast(z.data());}
  float* get_s() { return thrust::raw_pointer_cast(s.data());}
  float* get_p() { return thrust::raw_pointer_cast(p.data());}
  float* get_u() { return thrust::raw_pointer_cast(u.data());}
  float* get_x() { return thrust::raw_pointer_cast(x.data());}
  float* get_temp_res() { return thrust::raw_pointer_cast(temp_res.data());}
  void copyElementsFromHost(thrust::host_vector<float> v);
  void copyFlocalFromHost(thrust::host_vector<float> v);
private:
  thrust::device_vector<float> gpuB, gpuKlocals, gpuFlocals, gpuElements,
  diag, r, m, z, s, p, u, x, mask, n_adjelem, temp_res;
};

void TEST_THRUST();
void gpuAddWeighted(float *a, float *b, float v1, float v2, int size);
void gpuReductionWithMaskAndTransform(float *v, float *mask, int size, float *res, int size_new);
void gpuReductionWithMask(float *v, float *mask, int size, float *res, int size_new);
float gpuDotProduct(float *a, float *b, int size);


void gpuAddWeighted2(float *a_ptr, float *b_ptr,
                     float v1, float v2, int size);
void gpuTransformFrom_2N_to_6E(float *mask_ptr, float *a_ptr,
                               float *b_ptr, int size);
void gpuReductionWithMaskAndTransform2(float *v, float *mask, int size,
                                       float *res, int size_new);
void gpuReductionWithMask2(float *v, float *mask,
                           int size, float *res);
float gpuDotProduct2(float *a, float *b, int size);
void gpuDivideByElementwise(float *a, float *b,
                            float *c, int size);

void gpuCalculateKlocal2(gpuDataKeeper &gpu_data, int elementsCount,
                         float *h_nodesX, float *h_nodesY, int nodesCount,
                         float *h_D, float *h_constraints, int constraintsCount);

void gpuMultiplyKlocalByVec(gpuDataKeeper &gpu_data, int elementsCount);
void gpuCopyDeviceToDevice(float *data, float *dest, int size);
void gpuCopyDeviceToHost(float *data, float *dest, int size);
void gpuCountNAdjElem(gpuDataKeeper &gpu_data, int grid_size);

void gpuGenerateMask(gpuDataKeeper &gpuD, int elementsCount);
void copyElements(FEMdataKeeper FEMdata, gpuDataKeeper &gpuD);

#endif
