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

using namespace std;


void CalculateFEM(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE_vec(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateFEM_dyn(FEMdataKeeper &FEMdata, float rho, float alpha, float beta, float endtime, float dt, bool PRINT_DEBUG_INFO);
void CalculateFEM_dyn_vec(FEMdataKeeper &FEMdata, float rho, float damping_alpha, float damping_beta, float endtime, float dt, float beta1, float beta2, bool PRINT_DEBUG_INFO);

void gpuPCG_EbE_vec(FEMdataKeeper &FEMdata, MyArray &res, bool doAssemblyRes, float eps, bool PRINT_DEBUG_INFO);
void CalculateFEM_EbE_vec_GPU(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);
void CalculateNodeAdjElem(FEMdataKeeper FEMdata, std::unordered_map <int, std::vector<int>> &a);
void AssignLoadElement(FEMdataKeeper &FEMdata, std::unordered_map <int, std::vector<int>> nodeAdjElem);
void GetMapElement2Loadvector(FEMdataKeeper &FEMdata, std::unordered_map <int, MyArray> &loadVectors, float t);

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
                                   MyArray displacements);

SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata);
MyArray AssemblyF(FEMdataKeeper &FEMdata);


#endif
