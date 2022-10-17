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

using namespace std;


void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata);

void MakeResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper &FEMdata, ResultsDataKeeper &RESdata, std::string output_vtk, bool PRINT_DEBUG_INFO);
void FindConstraints(const std::vector<Constraint> constraints, std::vector<int> &indicesToConstraint);
void ApplyConstraints(SparseMatrixCOO& K, const std::vector<Constraint>& constraints, int n);

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix &D,
                                   std::vector<Element> &elements,
                                   MyArray &displacements,
                                   MyArray &all_B);

void CalculateFEM_EbE_vec_GPU(FEMdataKeeper &FEMdata, bool PRINT_DEBUG_INFO);

void ApplyConstraints_EbE(FEMdataKeeper &FEMdata);
//SparseMatrixCOO AssemblyStiffnessMatrix(FEMdataKeeper &FEMdata);


#endif
