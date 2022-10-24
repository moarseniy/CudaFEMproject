#ifndef FEMDATA_H
#define FEMDATA_H

#include <string.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

// https://stackoverflow.com/q/21667295
#include <regex>

#include "Linal2.h"
#include "Tools.h"
#include "femstruct.h"

using namespace std;

class FEMdataKeeper {
public:
  FEMdataKeeper(std::string name,
                int DIM,
                std::string project_directory,
                std::string prepared_meshes_directory,
                std::string results_directory,
                float poissonRatio, float youngModulus);

  void ShowInfo();
  std::string GetName()   { return this->name; }
  std::string GetResDir() { return this->res_dir; }

  int DIM;
  int nodesCount;
  int elementsCount;
  int boundaryEdgesCount;
  int loadsCount;
  int constraintsCount;

  Matrix D;
  std::vector<MyArray> nodes;

  MyArray CudaIndicesToConstraints;
  int CudaIndicesToConstraintsCount;

  std::vector<Element>        elements;
  std::vector<Constraint>     constraints;
  std::vector<Load>           loads;
  std::vector<BoundaryEdge>	  boundary;       // stores only boundary edges with nonzero b.c.
  MyArray displacements;
  MyArray pressure;

  MyArray all_B;

protected:
  void Parse();
  void ParseNodes();
  void ParseElements();
  void ParseConstraints();
  void ParseLoads();
  void ParseBoundaryEdges();
  void CreateMatrixD();

private:
  std::string name;
  std::string proj_dir;
  std::string prep_mesh_dir;
  std::string res_dir;
  float poissonRatio, youngModulus;
  fstream nodes_file, elements_file, loads_file, constraints_file, stress_file;

};

class ResultsDataKeeper {
public:
  ResultsDataKeeper(FEMdataKeeper *FEMdata,
                    bool withStressAlongAxis, bool withSmooth, bool withMises,
                    int nodesCount) {
    this->withStressAlongAxis = withStressAlongAxis;
    this->withMises = withMises;
    this->withSmooth = withSmooth;
    this->FEMdata = FEMdata;
    SmoothStress.Resize(nodesCount); //take it out somewhere else!!
  }
  //void MakeResults();

  bool withStressAlongAxis;
  bool withSmooth;
  bool withMises;

  std::vector<MyArray> Deformation;
  std::vector<MyArray> Stress;
  std::vector<float> sigma_mises;
  std::vector<float> epsilon_mises;
  std::vector<float> StressComponents;
  std::vector<float> MisesComponents;
  std::vector<float> StressComponentsSmooth;

  MyArray SmoothStress;
protected:
  FEMdataKeeper *FEMdata;

  //void CalculateStressAndDeformation();
};

#endif
