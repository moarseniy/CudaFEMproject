#ifndef FEMFUNC_H
#define FEMFUNC_H

#include "Linal2.h"
#include "Tools.h"

#include <iostream>
#include <vector>

using namespace std;


struct Element {
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
    void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B = Matrix(6, 12);
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
    Coords(int x, int y){
        this->x = x;
        this->y = y;
    }
    int x;
    int y;
};


class FEMdataKeeper {
public:
    FEMdataKeeper(){}

    void ParseFiles(std::string dir, std::string name, float poissonRatio, float youngModulus);
    void ShowInfo() {
        std::cout << "==========INFO==========" <<
                     "\nNodes count = " << nodesCount <<
                     "\nElements count = " << elementsCount <<
                     "\nConstraints count = " << constraintsCount <<
                     "\nLoads count = " << loadsCount <<
                     "\n========================\n\n";
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
    void SetDimensionCount(int dimension) {
        this->dimension = dimension;
    }

    void AllocateDynamicMemory() {
        nodesX.Resize(nodesCount);
        nodesY.Resize(nodesCount);
        nodesZ.Resize(nodesCount);
        loads.Resize(3 * nodesCount);
        D.Resize(6, 6);
    }

    int nodesCount;
    int elementsCount;
    int loadsCount;
    int constraintsCount;
    int dimension;

    Matrix D;
    MyArray nodesX;
    MyArray nodesY;
    MyArray nodesZ;
    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;
    MyArray loads;
};

void CalculateFiniteElementMethod(FEMdataKeeper &FEMdata, MyArray &displacements);


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
