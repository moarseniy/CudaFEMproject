#ifndef FEMFUNC_H
#define FEMFUNC_H

#include "Linal2.h"
#include <iostream>
#include <vector>

using namespace std;

struct Element
{
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
    void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B = Matrix(6, 12);
    int nodesIds[4];
};

struct ElementLight
{
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
    Coords(int x, int y){
        this->x = x;
        this->y = y;
    }
    int x;
    int y;
};

void FindConstraints(const std::vector<Constraint> constraints, std::vector<int> &indicesToConstraint);
void ApplyConstraints(SparseMatrixCOO& K, const std::vector<Constraint>& constraints, int n);

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<float> &epsilon_mises,
                                   std::vector<float> &sigma_mises,
                                   Matrix D,
                                   std::vector<Element> elements,
                                   MyArray displacements);


#endif
