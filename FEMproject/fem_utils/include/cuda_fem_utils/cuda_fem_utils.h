#ifndef CUDA_FEM_UTILS_H
#define CUDA_FEM_UTILS_H

#include <iostream>
#include <cassert>

//#include <cublas.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include <cuda_runtime.h>
#include <cuda.h>

#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/inner_product.h>

#include <fem_utils/fem_utils.h>
#include "cuda_fem_utils/fem_utils_kernels.h"


class CUDA_ElementsData : public ElementsData {
public:
  CUDA_ElementsData();
  CUDA_ElementsData(const dataKeeper &dk);
  ~CUDA_ElementsData();

  void getDiagonalElements() override;
  void transformWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) override;
  void applyConstraints() override;
  void calculateKlocal() override;
  void calculateKlocals() override;
  void calculateFlocal() override;
  void calculateFlocals() override;
  void calculateArea() override;
  void calculateLength() override;
  void genMask() override;
  void genAdjElements() override;
  void genCoordinates() override;
  void genFCoordinates() override;
  void genGradientMatrix() override;

};


#endif // CUDA_FEM_UTILS_H
