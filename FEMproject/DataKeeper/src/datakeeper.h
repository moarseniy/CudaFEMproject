#ifndef FEMDATA_H
#define FEMDATA_H

#include <string.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

// https://stackoverflow.com/q/21667295
#include <regex>

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

#include <Tools.h>
#include <femstruct.h>

#include <fstream>

struct MechanicalParams {
  // common params
  float poissonRatio, youngModulus, rho;
  // damping params
  float damping_alpha, damping_beta;
  // time params
  float dt, endtime;
  // newmark's scheme params
  float beta1, beta2;
};

struct dataPaths {
  // CudaFemProject path
  std::string proj_dir;
  // prepared_meshes directory path
  std::string prep_mesh_dir;

  // working directory path
  std::string work_dir;
  // results directory path
  std::string res_dir;
  // vtk file results path
  std::string output_vtk;

  std::string getResDir() { return this->res_dir; }
};


class dataKeeper {
public:
  dataKeeper(std::string configPath, std::string taskName);
  ~dataKeeper();

  void ShowInfo();
  std::string getTaskName() const { return taskName; }
  dataPaths getDataPaths() const { return paths; }
  size_t get_dim() const { return this->DIM; }
  size_t get_elementsCount() const { return elementsCount; }
  size_t get_nodesCount() const { return nodesCount; }
  size_t get_boundaryEdgesCount() const { return boundaryEdgesCount; }

  Matrix* get_nodes() const { return nodes; }
  Matrix* get_elementsIds() const { return elementsIds; }
  Matrix* get_constraintsIds() const { return constraintsIds; }
  Matrix* get_constraintsTypes() const { return constraintsTypes; }
  Matrix* get_boundaryAdjElems() const { return boundaryAdjElems; }
  Matrix* get_boundaryNodes() const { return boundaryNodes; }
  Matrix* get_boundaryNormals() const { return boundaryNormals; }
  Matrix* get_boundaryPressures() const { return boundaryPressureValues; }
  Matrix* get_D() const { return D; }
  Matrix* get_displacements() const { return displacements; }

protected:
  void parseJsonConfig(std::string jsonPath);
  void parseCounters();
  void allocateMemory();
  void CreateMatrixD();
  void ParseNodes();
  void ParseElements();
  void ParseConstraints();
  void ParseLoads();
  void ParseBoundaryEdges();
  void parseFiles();

  std::string taskName;
  std::string taskType;

  size_t DIM;
  size_t nodesCount;
  size_t elementsCount;
  size_t boundaryEdgesCount;
  size_t loadsCount;
  size_t constraintsCount;
  size_t constraintsIdsCounts;

  std::unordered_map <int, std::vector<int>> nodeAdjElem;

  Matrix *nodes;
  Matrix *elementsIds;
  Matrix *constraintsIds;
  Matrix *constraintsTypes;
  Matrix *D;

  Matrix *boundaryNormals;
  Matrix *boundaryAdjElems;
  Matrix *boundaryNodes;
  Matrix *boundaryPressureValues;

  std::vector<Element>        elements;
  std::vector<Constraint>     constraints;
  std::vector<Load>           loads;
  std::vector<BoundaryEdge>	  boundary;  // stores only boundary edges with nonzero b.c.

  Matrix *displacements;

  MechanicalParams mechParams;
  dataPaths paths;

  std::fstream nodes_file, elements_file, loads_file, constraints_file, stress_file;
};

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

  CPU_Matrix D;
  std::vector<CPU_Matrix> nodes;

  CPU_Matrix CudaIndicesToConstraints;
  int CudaIndicesToConstraintsCount;

  CPU_Matrix elementsNodesIds;
  CPU_Matrix all_B;

  std::vector<Element>        elements;
  std::vector<Constraint>     constraints;
  std::vector<Load>           loads;
  std::vector<BoundaryEdge>	  boundary;       // stores only boundary edges with nonzero b.c.
  Matrix *displacements;
  Matrix *pressure;

protected:
  void parseFiles();
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
  std::fstream nodes_file, elements_file, loads_file, constraints_file, stress_file;

};

class CUDA_FEMdataKeeper {

};

class ResultsDataKeeper {
public:
  ResultsDataKeeper(bool withStressAlongAxis, bool withSmooth, bool withMises,
                    int nodesCount) {
    this->withStressAlongAxis = withStressAlongAxis;
    this->withMises = withMises;
    this->withSmooth = withSmooth;
//    SmoothStress.Resize(nodesCount); //take it out somewhere else!!
  }
  //void MakeResults();

  bool withStressAlongAxis;
  bool withSmooth;
  bool withMises;

  std::vector<CPU_Matrix> Deformation;
  std::vector<CPU_Matrix> Stress;
  std::vector<float> sigma_mises;
  std::vector<float> epsilon_mises;
  std::vector<float> StressComponents;
  std::vector<float> MisesComponents;
  std::vector<float> StressComponentsSmooth;

  CPU_Matrix SmoothStress;

  //void CalculateStressAndDeformation();
};

#endif
