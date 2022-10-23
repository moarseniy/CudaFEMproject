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

void writeSnapshot(float t, int num_receivers, int grid_size, int n_gl_dofs, FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data);

void gpuCalculateFEM_DYN(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);

void gpuCalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN &gpu_data, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec_DYN_DAMP(FEMdataKeeper &FEMdata, gpuDataKeeper_DYN_DAMP &gpu_data, bool doAssemblyRes,  float eps, bool PRINT_DEBUG_INFO);

void CalculateNodeAdjElem(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> &a);
void AssignLoadElement(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem);
void GetMapElement2Loadvector(FEMdataKeeper &FEMdata, std::unordered_map <int, MyArray> &loadVectors, float t);
void ApplyConstraints_EbE(FEMdataKeeper &FEMdata);
void GenerateMask(FEMdataKeeper FEMdata, MyArray &mask);

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO);

void CalculateStressAlongAxis(std::vector<float> &StressComponents,
                              std::string axe,
                              std::string stress_component,
                              float fixed_value,
                              float a,
                              float b,
                              std::vector<MyArray> &Stress,
                              std::vector<MyArray> &nodes,
                              std::vector<Element> elements);

void CalculateStressAndDeformation(int DIM,
                                   std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix &D,
                                   std::vector<Element> &elements,
                                   MyArray &displacements, MyArray &all_B);

MyArray AssemblyF(FEMdataKeeper &FEMdata);


#endif
