#ifndef VTKfile_H
#define VTKfile_H

#include "Linal2.h"
#include "femfunc.h"
#include "Tools.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void MakeVTKfile2D(std::string output_vtk,
                   MyArray nodesX,
                   MyArray nodesY,
                   std::vector<Element> elements,
                   MyArray displacements,
                   std::vector<MyArray> Stress,
                   std::vector<float> sigma_mises,
                   std::vector<MyArray> Deformation,
                   std::vector<float> epsilon_mises,
                   MyArray SmoothStress);

void MakeVTKfile3D(std::string output_vtk,
                   MyArray nodesX,
                   MyArray nodesY,
                   MyArray nodesZ,
                   std::vector<Element> elements,
                   MyArray displacements,
                   std::vector<MyArray> Stress,
                   std::vector<float> sigma_mises,
                   std::vector<MyArray> Deformation,
                   std::vector<float> epsilon_mises);

#endif
