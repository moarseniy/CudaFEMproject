#ifndef FEMFUNC_H
#define FEMFUNC_H

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>
#include <Tools.h>
#include <datakeeper.h>
#include <femstruct.h>
#include <VTKfile.h>
#include <init.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>

void writeSnapshot(float t, int num_receivers, int grid_size, int n_gl_dofs, FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data);

void gpuCalculateFEM_DYN(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);

void gpuCalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, Matrix &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN_DAMP(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN_DAMP &gpu_data, bool doAssemblyRes,  float eps, bool PRINT_DEBUG_INFO);

void CalculateNodeAdjElem(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &a);
void AssignLoadElement(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem);
void GetMapElement2Loadvector(FEMdataKeeper &FEMdata, std::unordered_map <int, CPU_Matrix> &loadVectors, float t);
void ApplyConstraints_EbE(FEMdataKeeper &FEMdata);
void GenerateMask(FEMdataKeeper FEMdata, CPU_Matrix &mask);

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO);

void CalculateStressAlongAxis(std::vector<float> &StressComponents,
                              std::string axe,
                              std::string stress_component,
                              float fixed_value,
                              float a,
                              float b,
                              std::vector<CPU_Matrix> &Stress,
                              std::vector<CPU_Matrix> &nodes,
                              std::vector<Element> elements);

void CalculateStressAndDeformation(int DIM,
                                   std::vector<CPU_Matrix> &Deformation,
                                   std::vector<CPU_Matrix> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix &D,
                                   std::vector<Element> &elements,
                                   CPU_Matrix &displacements, CPU_Matrix &all_B);

CPU_Matrix AssemblyF(FEMdataKeeper &FEMdata);

// 2D only
void SmoothResults(std::string stress_component, CPU_Matrix &SmoothStress, std::vector<CPU_Matrix> Stress,
                   int nodesCount, std::vector<CPU_Matrix> &nodes, std::vector<Element> elements);
// 2D only
void CalculateStressAlongAxisSmooth(std::vector<float> &StressComponentsSmooth,
                                    std::string axe,
                                    float fixed_value,
                                    float a,
                                    float b,
                                    CPU_Matrix StressSmooth,
                                    std::vector<CPU_Matrix> &nodes,
                                    std::vector<Element> elements);
void CalculateMisesAlongLineMises(std::vector<float> &MisesComponents,
                                  float k, float m,
                                  float a, float b,
                                  std::vector<float> sigma_mises,
                                  std::vector<CPU_Matrix> &nodes,
                                  std::vector<Element> elements);

#endif
