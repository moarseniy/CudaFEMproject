#ifndef FEMFUNC_H
#define FEMFUNC_H

#include "Linal2.h"
#include "Tools.h"
#include "datakeeper.h"
#include "femstruct.h"
#include "VTKfile.h"
#include "init.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>

//// ------------------------------- CUDA --------------------------------------------
//#include <cuda_runtime.h>
//#include <cusparse_v2.h>

//#include <thrust/transform_reduce.h>
//#include <thrust/functional.h>
//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h>
//#include <thrust/inner_product.h>
//#include <thrust/device_ptr.h>

//#include <iostream>
//#include <cuda.h>
//#include <cusolverSp.h>
//// ---------------------------------------------------------------------------------


using namespace std;

void CalculateFEM(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_dyn(FEMdataKeeper &FEMdata, float rho, float alpha, float beta, float endtime, float dt, bool PRINT_DEBUG_INFO);
void CalculateFEM_dyn_vec(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);
void CalculateFEM_dyn_relaxation(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float dt, float beta1, float beta2, float eps);

void gpuCalculateFEM_DYN(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_implicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_implicit_DYN(FEMdataKeeper &FEMdata, float rho, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_implicit_DYN_DAMP(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_explicit(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_explicit_DYN(FEMdataKeeper &FEMdata, float rho, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_explicit_DYN_DAMP(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, bool PRINT_DEBUG_INFO);


void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN_DAMP(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN_DAMP &gpu_data, bool doAssemblyRes,  float eps, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE_vec_GPU(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);

void CalculateNodeAdjElem(FEMdataKeeper FEMdata, std::unordered_map <int, std::vector<int>> &a);
void AssignLoadElement(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem);
void GetMapElement2Loadvector(FEMdataKeeper &FEMdata, std::unordered_map <int, MyArray> &loadVectors, float t);
void solve_diag_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes);
void ApplyConstraints_EbE(FEMdataKeeper &FEMdata);
void GenerateMask(FEMdataKeeper FEMdata, MyArray &mask);

void MakeResults(FEMdataKeeper FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper FEMdata, ResultsDataKeeper RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO);

void FindConstraints(const std::vector<Constraint> constraints, std::vector<int> &indicesToConstraint);
void ApplyConstraints(SparseMatrixCOO& K, MyArray& F, const std::vector<Constraint>& constraints, int n);

void CalculateStressAlongAxis(std::vector<float> &StressComponents,
                              std::string axe,
                              std::string stress_component,
                              float fixed_value,
                              float a,
                              float b,
                              std::vector<MyArray> &Stress,
                              MyArray nodesX,
                              MyArray nodesY,
                              std::vector<Element> elements);

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix D,
                                   std::vector<Element> elements,
                                   MyArray displacements, MyArray all_B);

SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata);
MyArray AssemblyF(FEMdataKeeper &FEMdata);


#endif
