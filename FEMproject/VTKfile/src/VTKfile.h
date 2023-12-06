#ifndef VTKfile_H
#define VTKfile_H

#include <matrix_pack.h>

//#include <SegyIO.h>
//#include <SegyIn.h>
//#include <SegyOut.h>

#include <Tools.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void testSEGY(std::string path1, std::string path2);
void getSEGYinfo(std::string path);
void writeDisplForSEGY(std::string out_path, Matrix &src, std::vector<int> receivers, int axis, float timestep);
void convertToSEGY(std::string src_path, std::string SEGY_path, std::vector<int> receivers, float sampleInterval, int sampleNumber);

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
