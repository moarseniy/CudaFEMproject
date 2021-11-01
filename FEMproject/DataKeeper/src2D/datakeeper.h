#ifndef FEMDATA_H
#define FEMDATA_H

#include "Linal2.h"
#include "Tools.h"
#include "femstruct.h"
#include <iostream>
#include <vector>

using namespace std;



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

    void AllocateDynamicMemory() {
        nodesX.Resize(nodesCount);
        nodesY.Resize(nodesCount);
        loads.Resize(DIM * nodesCount);
        displacements.Resize(DIM * nodesCount);
        D.Resize(3 * (DIM - 1), 3 * (DIM - 1));
    }

    int nodesCount;
    int elementsCount;
    int loadsCount;
    int constraintsCount;

    Matrix D;
    MyArray nodesX;
    MyArray nodesY;
    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;
    MyArray loads;
    MyArray displacements;
};

#endif
