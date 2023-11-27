#ifndef FEMFUNC_H
#define FEMFUNC_H

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

#include <fem_utils/fem_utils.h>
#include <cuda_fem_utils/cuda_fem_utils.h>

#include <Tools.h>
#include <datakeeper.h>
#include <femstruct.h>
#include <VTKfile.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>

#include <segyUtils.h>


void gpuCalculateFEM_DYN(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);

void gpuCalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);

void gpuCalculateFEM_EbE_vec2(DEVICE_NAME deviceType, dataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec2(Matrix &displacements, DEVICE_NAME deviceType, dataKeeper &FEMdata, ElementsData &elemsData, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_DYN2(DEVICE_NAME devType, dataKeeper &dk, bool PRINT_DEBUG_INFO);

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, Matrix &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO);

void WriteResults(dataKeeper &dk, ResultsDataKeeper &RESdata);

#endif
