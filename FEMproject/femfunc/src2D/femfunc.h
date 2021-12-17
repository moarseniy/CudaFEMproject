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

using namespace std;


void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata);

void MakeResults(FEMdataKeeper &FEMdata, std::string output_vtk);

void FindConstraints(const std::vector<Constraint> constraints, std::vector<int> &indicesToConstraint);
void ApplyConstraints(SparseMatrixCOO& K, MyArray& F, const std::vector<Constraint>& constraints, int n);

void CalculateStressAlongAxe(std::vector<float> &StressComponents,
                             std::string axe,
                             std::string stress_component,
                             float fixed_value,
                             float a,
                             float b,
                             std::vector<MyArray> Stress,
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
