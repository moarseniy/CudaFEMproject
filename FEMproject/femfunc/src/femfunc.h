#ifndef FEMFUNC_H
#define FEMFUNC_H

#include <matrix_pack.h>
#include <fem_utils.h>

#include <Tools.h>
#include <datakeeper.h>
#include <femstruct.h>
#include <VTKfile.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>

//#include <segyUtils.h>

void gpuCalculateFEM_EbE_vec(DEVICE_NAME deviceType, dataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void gpuPCG_EbE_vec(Matrix &displacements, DEVICE_NAME deviceType, dataKeeper &FEMdata, ElementsData &elemsData, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void gpuCalculateFEM_DYN(DEVICE_NAME devType, dataKeeper &dk, bool PRINT_DEBUG_INFO);

void WriteResults(dataKeeper &dk, ResultsDataKeeper &RESdata);

#endif
