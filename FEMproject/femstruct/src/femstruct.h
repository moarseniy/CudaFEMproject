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

using namespace std;

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
  int node0, node1, node2;
  int adj_elem1, adj_elem2;   // if there is no second adjacement element, set adj_elem2 equaled to -1
  // see BoundaryEdge structure
};

struct BoundaryEdge : Edge, TimeDependentEntity {
  BoundaryEdge() {
    adj_elem2 = -1;
    type0 = Constraint::Type::NO;
    type1 = Constraint::Type::NO;
    type2 = Constraint::Type::NO;
  }
  float normal_x, normal_y, normal_z;   // make sure it is normalized
  Constraint::Type type0, type1, type2;
};

struct Element {
  // CalculateStiffnessMatrix will be deprecated!
  void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY);
  void CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
  void CalculateMlocal(float rho, MyArray& nodesX, MyArray& nodesY, bool lumped);
  void CalculateClocal(float alpha, float beta);
  void CalculateFlocal(BoundaryEdge& edge, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ, float t);
  int Global2LocalNode(int glob_num);
  Element() {
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
    m.Resize(6 * (DIM - 1)); x.Resize(6 * (DIM - 1)); r.Resize(6 * (DIM - 1));
    z.Resize(6 * (DIM - 1)); s.Resize(6 * (DIM - 1)); p.Resize(6 * (DIM - 1));
    u.Resize(6 * (DIM - 1));

    // for dynamics
    x_pred.Resize(6 * (DIM - 1)); vel_pred.Resize(6 * (DIM - 1)); vel.Resize(6 * (DIM - 1));

    // for debug in PCG EbE. Delete after is done!
    b.Resize(6 * (DIM - 1));
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

  int nodesIds[4];
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

  void assignElement(std::unordered_map <int, std::vector<int>> nodeAdjElem);
};

struct Couple {
  int index;
  int degree;
};




#endif
