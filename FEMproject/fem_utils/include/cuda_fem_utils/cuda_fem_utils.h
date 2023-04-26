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

  void getDiagonalElements(Matrix &Locals, Matrix &tgt) override;
  void transformWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMask(Matrix &src, Matrix &dest) override;
  void reductionWithMaskAndTransform(Matrix &src, Matrix &dest, size_t size) override;
  void applyConstraints() override;
  void calculateKlocal() override;
  void calculateKlocals() override;
  void calculateFlocal(float t, const WaveletParams &waveParams) override;
  void calculateFlocals(float t, const WaveletParams &waveParams) override;
  void calculateArea() override;
  void calculateLength() override;
  void calculateLength3D() override;
  void genMask() override;
  void genAdjElements() override;
  void genCoordinates() override;
  void genFCoordinates() override;
  void genGradientMatrix() override;
  void genGradientMatrix3D() override;

  void calculateDiag(Matrix &diag, float cM = 0.f, float cK = 0.f, float cC = 0.f, float dampAlpha = 0.f, float dampBeta = 0.f) override;
  void solveDiagSystem(Matrix &diagonal, Matrix &v, Matrix &tgt, bool transformRes) override;
  void calculateMlocals(bool isLumped, const MechanicalParams &mechParams) override;
};


#endif // CUDA_FEM_UTILS_H
