#ifndef FEMSTRUCT_H
#define FEMSTRUCT_H

#include "Linal2.h"
#include "Tools.h"
//#include "femfunc.h"

#include <iostream>
#include <vector>

using namespace std;

struct Edge {
    int node0, node1;
    int adj_elem1, adj_elem2;   // if there is no second adjacement element, set adj_elem2 equaled to -1
                                // see BoundaryEdge structure
};

struct BoundaryEdge : Edge {
    BoundaryEdge() {
        adj_elem2 = -1;
    }
    float normal_x, normal_y;   // make sure it is normalized
};

struct Element {
    // CalculateStiffnessMatrix will be deprecated!
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY);
    void CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY);
    void CalculateFlocal(BoundaryEdge& edge, MyArray& nodesX, MyArray& nodesY,  float pressure_value);
    Element() {
        B.Resize(3 * (DIM - 1), 6 * (DIM - 1));
        Klocal.Resize(6 * (DIM - 1), 6 * (DIM - 1));
        Flocal.Resize(6 * (DIM - 1));
        //Q.Resize();
    }
    //void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B;
    Matrix Klocal;
    MyArray Flocal;
    //Matrix Q;
    int nodesIds[3];
    //float jacobian; // node 1 -> (0,0) ; node 2 -> (1,0) ; node 3 -> (0,1)
};

struct ElementLight {
    int nodesIds[3];
};

struct Constraint {
    enum Type {
        UX = 1,
        UY = 2,
        UZ = 3
    };
    int node;
    Type type;
};

struct Couple {
    int index;
    int degree;
};

struct Coords {
    Coords(int x, int y) {
        this->x = x;
        this->y = y;
    }
    int x;
    int y;
};

typedef std::vector<Triplet> vecSparseMatrixCOO;


#endif
