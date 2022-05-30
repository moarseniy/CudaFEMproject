#ifndef FEMDATA_H
#define FEMDATA_H

#include <string.h>
#include <iostream>
#include <vector>
#include <sstream>

// https://stackoverflow.com/q/21667295
#include <regex>

#include "Linal2.h"
#include "Tools.h"
#include "femstruct.h"

using namespace std;

class FEMdataKeeper {
public:
  FEMdataKeeper(std::string name,
                std::string project_directory,
                std::string prepared_meshes_directory,
                std::string results_directory) {
    this->name = name;
    this->proj_dir = project_directory;
    this->prep_mesh_dir = prepared_meshes_directory;
    this->res_dir = results_directory;

    this->nonzeroMatrixNumCount = -1;
  }

  void ParseFiles(float poissonRatio, float youngModulus);
  void ShowInfo() {
    std::cout << "==========INFO==========" <<
                 "\nTask name = " << name <<
                 "\nNodes count = " << nodesCount <<
                 "\nElements count = " << elementsCount <<
                 "\nConstraints count = " << constraintsCount <<
                 "\nLoads count = " << loadsCount <<
                 "\nBoundary edges with b.c. count = " << boundaryEdgesCount <<
                 "\n========================\n\n";
  }
  std::string get_name() {
    return this->name;
  }

  void SetNodesCount(int nodesCount) {
    this->nodesCount = nodesCount;
  }
  void SetElementsCount(int elementsCount) {
    this->elementsCount = elementsCount;
  }
  void SetLoadsCount(int loadsCount) {
    this->loadsCount = loadsCount;
  }
  void SetConstraintsCount(int constraintsCount) {
    this->constraintsCount = constraintsCount;
  }
  void SetBoundaryEdgesCount(int boundaryEdgesCount) {
    this->boundaryEdgesCount = boundaryEdgesCount;
  }

  void AllocateDynamicMemory() {
    nodesX.Resize(nodesCount);
    nodesY.Resize(nodesCount);
    //loads.Resize(DIM * nodesCount);
    displacements.Resize(DIM * nodesCount);
    D.Resize(3 * (DIM - 1), 3 * (DIM - 1));
    pressure.Resize(boundaryEdgesCount);
    all_B.Resize(3 * 6 * elementsCount);
  }

  std::string name;
  std::string proj_dir;
  std::string prep_mesh_dir;
  std::string res_dir;

  int nodesCount;
  int elementsCount;
  int boundaryEdgesCount;
  int loadsCount;
  int constraintsCount;
  int nonzeroMatrixNumCount;

  Matrix D;
  MyArray nodesX;
  MyArray nodesY;

  MyArray CudaIndicesToConstraints;
  int CudaIndicesToConstraintsCount;

  std::vector<Element>        elements;
  std::vector<Constraint>     constraints;
  std::vector<Load>           loads;
  std::vector<BoundaryEdge>	  boundary;       // stores only those boundary edges
  // that have nonzero boundary conditions
  //MyArray loads;
  MyArray displacements;
  MyArray pressure;

  MyArray all_B;

};

class ResultsDataKeeper {
public:
  ResultsDataKeeper(bool withSmooth, bool withMises, int nodesCount) {
    this->withMises = withMises;
    this->withSmooth = withSmooth;
    SmoothStress.Resize(nodesCount); //take it out somewhere else!!
  }

  void AllocateDynamicMemory() {
    //SmoothStress.Resize();
  }

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
};

#endif
