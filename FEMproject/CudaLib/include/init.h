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
  float* get_tmp() { return thrust::raw_pointer_cast(tmp.data());}
  float* get_loads() { return thrust::raw_pointer_cast(loads.data());}
  float* get_temp_res() { return thrust::raw_pointer_cast(temp_res.data());}

  void copyElementsFromHost(thrust::host_vector<float> v);
  void copyFlocalFromHost(thrust::host_vector<float> v);
  void copyLoadsFromHost(thrust::host_vector<float> v);
  void copyBmatrixToHost(float *all_B);
  void setZeroVec();
protected:
  thrust::device_vector<float> gpuB, gpuKlocals, gpuFlocals, gpuElements,
  diag, r, m, z, s, p, u, x, mask, n_adjelem, temp_res, loads;
  thrust::device_vector<float> tmp; // for any temporarily calculations if needed
};

struct WeightedAddCoef {
//  WeightedAddCoef(float val1, float val2);
//  WeightedAddCoef(float val1, float val2, float val3);
  float val1, val2, val3;
};

class gpuDataKeeper_DYN : public gpuDataKeeper {
public:
  gpuDataKeeper_DYN(int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped);
  ~gpuDataKeeper_DYN();
  float* get_Mlocals() { return thrust::raw_pointer_cast(gpuMlocals.data());}
  float* get_diagM() { return thrust::raw_pointer_cast(diagM.data());}
  float* get_displ() { return thrust::raw_pointer_cast(displ.data());}
  float* get_displ_global() { return thrust::raw_pointer_cast(displ_global.data());}
  float* get_vel() { return thrust::raw_pointer_cast(vel.data());}


  void set_SLAU_matrix_coefs(float cM, float cK);

  WeightedAddCoef SLAU_matrix_coefs; // A = val1*M + val2*K + val3*C
  bool isLumped;
protected:
  // acceleration (res) will be stored in x
  thrust::device_vector<float> gpuMlocals, vel, diagM, displ, displ_global;
};

class gpuDataKeeper_DYN_DAMP : public gpuDataKeeper_DYN {
public:
  gpuDataKeeper_DYN_DAMP(int elementsCount, int nodesCount, bool doAssemblyRes, bool isLumped, float damping_alpha, float damping_beta);
  ~gpuDataKeeper_DYN_DAMP();
  void set_SLAU_matrix_coefs(float cM, float cK, float cC);
  void set_damping_coefs(float cM, float cK);

  WeightedAddCoef damping_coefs; // C = val1*M + val2*K
};

void TEST_THRUST();
void TEST_THRUST_GRISHA();
void TEST_THRUST_GenerateMask(FEMdataKeeper &FEMdata);
void TEST_THRUST_MultiplyByVector();
void TEST_THRUST_thrustReductionWithMask();
void TEST_THRUST_thrustTransform_2N_to_6E();
void TEST_THRUST_CNorm();
void gpuAddWeighted(float *a, float *b, float v1, float v2, int size);
void gpuReductionWithMaskAndTransform(float *v, float *mask, int size, float *res, int size_new);
void gpuReductionWithMask(float *v, float *mask, int size, float *res, int size_new);
float gpuDotProduct(float *a, float *b, int size);
void gpuCalculateFEM_dyn_explicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1);
void gpuCalculateFEM_dyn_relaxation(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float dt, float beta1, float eps);


void gpuAddWeighted2(float *a_ptr, float *b_ptr,
                     float v1, float v2, int size);
void gpuAdd(float *a_ptr, float *b_ptr, int size);
void gpuDivide(float *a_ptr, float *b_ptr, int size);
void gpuDivide_res(float *a_ptr, float *b_ptr, int size, float *res);
void gpuCopy(float *v, float *dest, int size);
void gpuTransformFrom_2N_to_6E(float *mask_ptr, float *a_ptr,
                               float *b_ptr, int size);
void gpuReductionWithMaskAndTransform2(float *v, float *mask, int size,
                                       float *res, int size_new);
void gpuReductionWithMask2(float *v, float *mask,
                           int size, float *res);
float gpuDotProduct2(float *a, float *b, int size);
void gpuDivideByElementwise(float *a, float *b,
                            float *c, int size);

void gpuCalculateKlocal2(gpuDataKeeper &gpu_data, FEMdataKeeper &FEMdata);
void gpuCalculateMlocal(gpuDataKeeper &gpu_data, int elementsCount,
                         float *h_nodesX, float *h_nodesY, int nodesCount);

void gpuMultiplyKlocalByVec(gpuDataKeeper &gpu_data, int elementsCount);
void gpuMultiplyAlocalByVec(gpuDataKeeper_DYN &gpu_data, int elementsCount);
void gpuMultiplyAlocalByVec_DAMP(gpuDataKeeper_DYN_DAMP &gpu_data, int elementsCount);
void gpuMultiplyMatrixByVec(float* Matr, float* Vec, float* Res, int elementsCount);
void gpuMultiplyClocalByVec(gpuDataKeeper_DYN_DAMP &gpu_data, float* Vec, float* Res, int elementsCount);
void gpuCalculateMlocal(gpuDataKeeper_DYN &gpu_data, int elementsCount,
                         float *h_nodesX, float *h_nodesY, int nodesCount, float rho);
void gpuCalculateDiag(gpuDataKeeper_DYN &gpu_data, int elementsCount);
void gpuCalculateDiag_DAMP(gpuDataKeeper_DYN_DAMP &gpu_data, int elementsCount);
void gpuCopyDeviceToDevice(float *data, float *dest, int size);
void gpuCopyDeviceToHost(float *data, float *dest, int size);
void gpuCopyHostToDevice(float *data, float *dest, int size);
void gpuCountNAdjElem(gpuDataKeeper &gpu_data, int grid_size);

void gpuGenerateMask(gpuDataKeeper &gpuD, int elementsCount);
void gpuTransform_2N_to_6E(float *d_v, int n_gl_dofs, float *d_mask, float *d_res, int grid_size);
void gpuSolveDiag(float *diag, float *r, float *res,
                  float *mask, float *n_adjelem,
                  int grid_size, int n_gl_dofs,
                  bool doAssemblyRes);
float gpuCNorm(float *v, int size);
void copyElementsAndFlocals(FEMdataKeeper &FEMdata, gpuDataKeeper &gpuD);
void copyFlocals(FEMdataKeeper &FEMdata, gpuDataKeeper &gpuD);
void copyLoads(gpuDataKeeper &gpuD, std::unordered_map <int, MyArray> &loadVectors, int elementsCount);

#endif
