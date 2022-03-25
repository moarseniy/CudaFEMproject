#ifndef FEMFUNC_H
#define FEMFUNC_H

#include "Linal2.h"
#include "Tools.h"
#include "datakeeper.h"
#include "femstruct.h"
#include "VTKfile.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <algorithm>

using namespace std;


void CalculateFEM(FEMdataKeeper &FEMdata);
void CalculateFEM_EbE(FEMdataKeeper &FEMdata);
void CalculateNodeAdjElem(FEMdataKeeper FEMdata, std::unordered_map <int, std::vector<int>> &a);

void MakeResults(FEMdataKeeper FEMdata, ResultsDataKeeper &RESdata);
void WriteResults(FEMdataKeeper FEMdata, ResultsDataKeeper RESdata, std::string output_vtk);

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
