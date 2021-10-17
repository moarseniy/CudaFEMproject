#ifndef FEMSTRUCT_H
#define FEMSTRUCT_H

#include "Linal2.h"
#include "Tools.h"
//#include "femfunc.h"

#include <iostream>
#include <vector>

using namespace std;

struct Element {
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
    Element() {
        B.Resize(3 * (DIM - 1), 6 * (DIM - 1));
    }
    //void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B;
    int nodesIds[4];
};

struct ElementLight {
    int nodesIds[4];
};

struct Constraint {
    enum Type {
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


#endif
