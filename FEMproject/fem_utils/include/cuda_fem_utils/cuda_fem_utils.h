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

  void getDiagonalElements(Matrix &K, size_t el) override;
  void transformWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) override;
  void applyConstraints(size_t el) override;
  void calculateKlocal(size_t elementId) override;
  void calculateKlocals() override;
  size_t getLocalId(size_t elementId, size_t nodeId) override;
  void calculateFlocal(size_t elementId) override;
  void calculateFlocals() override;
  float calculateArea(size_t elementId) override;
  float calculateLength(size_t elementId) override;
  void genMask() override;
  void genAdjElements() override;
  void genCoordinates(size_t elementId) override;
  void genFCoordinates(size_t elementId) override;
  float genGradientMatrix(size_t elementId) override;

};


#endif // CUDA_FEM_UTILS_H
