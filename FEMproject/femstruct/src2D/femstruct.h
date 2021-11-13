#ifndef FEMSTRUCT_H
#define FEMSTRUCT_H

#include "Linal2.h"
#include "Tools.h"
//#include "femfunc.h"

#include <iostream>
#include <vector>

using namespace std;

struct Element {
    // CalculateStiffnessMatrix will be deprecated!
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY);
    void CalculateKlocal(Matrix& D, MyArray& nodesX, MyArray& nodesY);
    Element() {
        B.Resize(3 * (DIM - 1), 6 * (DIM - 1));
        Klocal.Resize(6 * (DIM - 1), 6 * (DIM - 1));
        //Q.Resize();
    }
    //void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B;
    Matrix Klocal;
    //Matrix Q;
    int nodesIds[3];
    float jacobian; // node 1 -> (0,0) ; node 2 -> (1,0) ; node 3 -> (0,1)
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
