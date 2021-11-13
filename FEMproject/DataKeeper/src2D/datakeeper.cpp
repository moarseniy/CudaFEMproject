
#include "datakeeper.h"


using namespace std;

void FEMdataKeeper::ParseFiles(std::string dir, std::string name, float poissonRatio, float youngModulus) {
    CheckRunTime(__func__)

    fstream nodes_file, elements_file, loads_file, constraints_file;

    nodes_file.open(dir + name + "/nodes.txt", fstream::in);
    elements_file.open(dir + name + "/elements.txt", fstream::in);
    loads_file.open(dir + name + "/loads.txt", fstream::in);
    constraints_file.open(dir + name + "/constraints.txt", fstream::in);

    nodes_file >> nodesCount;
    elements_file >> elementsCount;
    constraints_file >> constraintsCount;
    loads_file >> loadsCount;

    AllocateDynamicMemory();

    D(0,0) = 1.0;			D(0, 1) = poissonRatio;	D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0; 			D(1, 2) = 0.0;
    D(2, 0) = 0.0;        	D(2, 1) = 0.0;        	D(2, 2) = (1.0 - poissonRatio) / 2.0;
    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));

    for (int i = 0; i < nodesCount; ++i) {
        nodes_file >> nodesX[i] >> nodesY[i];
    }

    for (int i = 0; i < elementsCount; ++i) {
        Element element;
        elements_file >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
        float X2 = nodesX[element.nodesIds[1]] - nodesX[element.nodesIds[0]];
        float X3 = nodesX[element.nodesIds[2]] - nodesX[element.nodesIds[0]];
        float Y2 = nodesY[element.nodesIds[1]] - nodesY[element.nodesIds[0]];
        float Y3 = nodesY[element.nodesIds[2]] - nodesY[element.nodesIds[0]];
        element.jacobian = 1.0 / (X3 * Y2 - X2 * Y3);
        elements.push_back(element);
    }

    for (int i = 0; i < constraintsCount; ++i) {
        Constraint constraint;
        int type;
        constraints_file >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints.push_back(constraint);
    }

    for (int i = 0; i < loadsCount; ++i) {
        int node;
        float x, y;
        loads_file >> node >> x >> y;
        loads[2 * node + 0] = x;
        loads[2 * node + 1] = y;
    }
}
