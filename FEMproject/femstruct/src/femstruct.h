#ifndef FEMSTRUCT_H
#define FEMSTRUCT_H

#include <iostream>
#include <vector>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <math.h>

#include "Linal2.h"
#include "Tools.h"
//#include "femfunc.h"

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

struct Element {
  // CalculateStiffnessMatrix will be deprecated!
  void CalculateKlocal(Matrix& D, std::vector<MyArray> &nodes);
  void CalculateMlocal(float rho, std::vector<MyArray> &nodes, bool lumped);
  void CalculateClocal(float alpha, float beta);
  void CalculateFlocal2D(BoundaryEdge& edge, std::vector<MyArray> &nodes, float t);
  void CalculateFlocal3D(BoundaryEdge& edge, std::vector<MyArray> &nodes, float t);
  int Global2LocalNode(int glob_num);
  Element(int DIM) {
    this->DIM = DIM;
    nodesIds.resize(DIM + 1);
    B.Resize(3 * (DIM - 1), 6 * (DIM - 1));
    Klocal.Resize(6 * (DIM - 1), 6 * (DIM - 1));
    Mlocal.Resize(6 * (DIM - 1), 6 * (DIM - 1));
    Clocal.Resize(6 * (DIM - 1), 6 * (DIM - 1)); alpha = 0; beta = 0;
    Alocal.Resize(6 * (DIM - 1), 6 * (DIM - 1));
    Flocal.Resize(6 * (DIM - 1));
    blocal.Resize(6 * (DIM - 1));
    res.Resize(6 * (DIM - 1));
    //Q.Resize();

    // for PCG EbE
    m.Resize(3 * DIM); x.Resize(3 * DIM); r.Resize(3 * DIM);
    z.Resize(3 * DIM); s.Resize(3 * DIM); p.Resize(3 * DIM);
    u.Resize(3 * DIM);

    // for dynamics
    x_pred.Resize(3 * DIM); vel_pred.Resize(3 * DIM); vel.Resize(3 * DIM);

    // for debug in PCG EbE. Delete after is done!
    b.Resize(3 * DIM);
  }
  //void FindSparseSize(std::vector<couple> &Sparse);
  Matrix B;
  Matrix Klocal;
  Matrix Mlocal;
  Matrix Clocal; float alpha, beta;   // C = alpha * M + beta * K;
  MyArray Flocal;

  Matrix Alocal;
  MyArray blocal;
  MyArray res;
  //Matrix Q;

  // for dynamics
  // --------------------------
  MyArray x_pred, vel_pred, vel;
  // --------------------------

  // for PCG EbE
  // --------------------------
  MyArray m, x, r, z, s, p, u;
  // --------------------------

  // for debug in PCG EbE. Delete after is done!
  // --------------------------
  MyArray b;
  // --------------------------

  int DIM;
  std::vector<int> nodesIds;

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
