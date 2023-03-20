#ifndef FEMSTRUCT_H
#define FEMSTRUCT_H

#include <iostream>
#include <vector>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <math.h>

#include <matrix_pack/matrix_pack.h>
#include <cuda_matrix_pack/cuda_matrix.h>

#include <Tools.h>

// https://stackoverflow.com/a/313990
#include <algorithm>
#include <cctype>
#include <string>

struct TimeDependentEntity {
  float value;
  float ampl, freq, timeshift;

  void (TimeDependentEntity::*wavelet)(float);

  TimeDependentEntity() {
    wavelet = &TimeDependentEntity::Constant;
    timeshift = 0.0f;
    ampl = 0.0f;    // no stress applied
  }

  bool isFloat(const std::string& str);
  void parseString(std::string& str);
  void update(float t);
  void Ricker(float t);
  void Berlage(float t);
  void Constant(float t);
};

struct Constraint {
  enum Type {
    NO = 0,
    UX = 1,
    UY = 2,
    UZ = 3,
    UXY = 4,
    UXZ = 5,
    UYZ = 6,
    UXYZ = 7
  };
  int node;
  Type type;
};

struct Edge {
  Edge(){}
  Edge(int DIM) {
    node.resize(DIM);
  }
  std::vector<int> node;
  int adj_elem1, adj_elem2;   // if there is no second adjacement element, set adj_elem2 equaled to -1
  // see BoundaryEdge structure
};

struct BoundaryEdge : Edge, TimeDependentEntity {
  BoundaryEdge(int DIM) {
    adj_elem2 = -1;
    for (int i = 0; i < DIM; ++i) {
      type.push_back(Constraint::Type::NO);
    }
    normal.resize(DIM);
    node.resize(DIM);
  }
  std::vector<float> normal;   // make sure it is normalized
  std::vector<Constraint::Type> type;
};

//class ElementsData {
//public:
//  ElementsData(size_t DIM, size_t elementsCount);
//  void CalculateKlocal();

//private:
//  Matrix *Flocals;
//  Matrix *Klocals;
//  Matrix *B;
//};

struct Element {
  float calculateArea(CPU_Matrix &Coordinates);
  void calculateGradientMatrix(CPU_Matrix &Coordinates, float area);
  void calculateCoordinates(CPU_Matrix &Coordinates,
                            std::vector<Matrix> &nodes);

  // CalculateStiffnessMatrix will be deprecated!
  void CalculateKlocal(Matrix& D, std::vector<CPU_Matrix> &nodes);
  void CalculateMlocal(float rho, std::vector<CPU_Matrix> &nodes, bool lumped);
  void CalculateClocal(float alpha, float beta);
  void CalculateFlocal2D(BoundaryEdge& edge, std::vector<CPU_Matrix> &nodes, float t);
  void CalculateFlocal3D(BoundaryEdge& edge, std::vector<CPU_Matrix> &nodes, float t);
  int Global2LocalNode(int glob_num);
  Element(int dim);
  ~Element();

  int DIM;
  std::vector<int> nodesIds;

  CPU_Matrix B;
  CPU_Matrix Klocal;
  CPU_Matrix Mlocal;
  CPU_Matrix Clocal;
  float alpha, beta;   // C = alpha * M + beta * K;
  CPU_Matrix Flocal;

  CPU_Matrix Alocal;
  CPU_Matrix blocal;
  CPU_Matrix res;

  // for dynamics
  // --------------------------
  CPU_Matrix x_pred, vel_pred, vel;
  // --------------------------

  // for PCG EbE
  // --------------------------
  CPU_Matrix m, x, r, z, s, p, u;
  // --------------------------

  // for debug in PCG EbE. Delete after is done!
  // --------------------------
  CPU_Matrix b;
  // --------------------------

  //float jacobian; // node 1 -> (0,0) ; node 2 -> (1,0) ; node 3 -> (0,1)
};

struct Load : TimeDependentEntity {
  int dof;
  int elem;   // For EbE. If the node has two or more adjacent element,
  // then only one of these elements has the total load contribution

  Load() {
    TimeDependentEntity();
    elem = -1; // no element assigned
  }

  void assignElement(int DIM, std::unordered_map <int, std::vector<int>> &nodeAdjElem);
};

#endif
