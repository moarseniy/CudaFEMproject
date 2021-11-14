#ifndef VTKfile_H
#define VTKfile_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Linal2.h"
#include "femfunc.h"
#include "Tools.h"

void MakeVTKfile(std::string output_vtk,
                 MyArray nodesX,
                 MyArray nodesY,
                 std::vector<Element> elements,
                 MyArray displacements,
                 std::vector<MyArray> Stress,
                 std::vector<float> sigma_mises,
                 std::vector<MyArray> Deformation,
                 std::vector<float> epsilon_mises);

void MakeVTKfileNew(std::string output_vtk,
                 MyArray nodesX,
                 MyArray nodesY,
                 std::vector<ElementLight> elements,
                 MyArray displacements);

#endif
