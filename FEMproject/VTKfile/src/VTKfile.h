#ifndef VTKfile_H
#define VTKfile_H

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

#include "femfunc.h"
#include "Tools.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void MakeVTKfile2D(std::string output_vtk,
                   std::vector<CPU_Matrix> &nodes,
                   std::vector<Element> &elements,
                   Matrix &displacements,
                   std::vector<CPU_Matrix> &Stress,
                   std::vector<float> &sigma_mises,
                   std::vector<CPU_Matrix> &Deformation,
                   std::vector<float> &epsilon_mises,
                   CPU_Matrix &SmoothStress);

void MakeVTKfile3D(std::string output_vtk,
                   std::vector<CPU_Matrix> &nodes,
                   std::vector<Element> &elements,
                   Matrix &displacements,
                   std::vector<CPU_Matrix> &Stress,
                   std::vector<float> &sigma_mises,
                   std::vector<CPU_Matrix> &Deformation,
                   std::vector<float> &epsilon_mises);

void MakeVTKfile2D(std::string output_vtk,
                   Matrix &nodes,
                   Matrix &elements,
                   Matrix &displacements);

void MakeVTKfile3D(std::string output_vtk,
                   Matrix &nodes,
                   Matrix &elements,
                   Matrix &displacements);

#endif
