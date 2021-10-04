#include "Linal2.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
//#include <cmath>

#include "femfunc.h"
#include "VTKfile.h"
#include "Parser.h"

using namespace std;

int main(void) {

    std::string name = "bulk_test";


    std::string directory = "C:/Users/mokin/Desktop/FiniteElementMethod3D/prepared_meshes/";
    std::string output_vtk = "C:/Users/mokin/Desktop/FiniteElementMethod3D/final_results/results.vtk";
    std::string output_results = "C:/Users/mokin/Desktop/FiniteElementMethod3D/final_results/output.txt";

    fstream nodes_file, elements_file, loads_file, constraints_file;
    nodes_file.open(directory + name + "/nodes.txt", fstream::in);
    elements_file.open(directory + name + "/elements.txt", fstream::in);
    loads_file.open(directory + name + "/loads.txt", fstream::in);
    constraints_file.open(directory + name + "/constraints.txt", fstream::in);


    int start_time = clock();
    unsigned long long int end_time;

    float poissonRatio = 0.3, youngModulus = 2e+11;
    Matrix D(6, 6);

    D(0, 0) = D(1, 1) = D(2, 2) = 1.0;
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = poissonRatio / (1.0 - poissonRatio);
    D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));

    D.scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));

    int nodesCount;
    nodes_file >> nodesCount;

    MyArray nodesX(nodesCount);
    MyArray nodesY(nodesCount);
    MyArray nodesZ(nodesCount);

    for (int i = 0; i < nodesCount; ++i) {
        nodes_file >> nodesX[i] >> nodesY[i] >> nodesZ[i];
    }

    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;

    int elementCount;
    elements_file >> elementCount;
    for (int i = 0; i < elementCount; ++i) {
        Element element;
        elements_file >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2] >> element.nodesIds[3];
        elements.push_back(element);
    }

    int constraintCount;
    constraints_file >> constraintCount;

    for (int i = 0; i < constraintCount; ++i) {
        Constraint constraint;
        int type;
        constraints_file >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints.push_back(constraint);
    }

    MyArray loads(3 * nodesCount);

    int loadsCount;
    loads_file >> loadsCount;

    for (int i = 0; i < loadsCount; ++i) {
        int node;
        float x, y, z;
        loads_file >> node >> x >> y >> z;
        loads[3 * node + 0] = x;
        loads[3 * node + 1] = y;
        loads[3 * node + 2] = z;
    }
    //loads.Show();

    cout << "Data read success\n";  

    end_time = clock();
    cout<< "Time: "<< end_time - start_time<< " ms\n";

    ///////////////////////////

    std::vector<Triplet> triplets;
    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        it->CalculateStiffnessMatrix(D, triplets, nodesX, nodesY, nodesZ);
    }


    cout << "\nCalculateStiffnessMatrix success\n";

    end_time = clock();
    cout << "Time(CalculateStiffMat): "<< end_time - start_time<< " ms\n";

    cout << "Triplets Size = " << triplets.size() << endl;

    MyArray displacements(loads.get_size());
    SparseMatrixCOO globalK(triplets.size());

    globalK.ConvertTripletToSparse(triplets);
    cout << "new size= "<<globalK.get_size()<<"\n";
    cout << "!!!" << sizeof (globalK) << endl;

    end_time = clock();
    cout<< "Time(ConvertTripletToSparse): "<< end_time - start_time<< " ms\n";

    globalK.resize();

    end_time = clock();
    cout<< "Time(resize): "<< end_time - start_time<< " ms\n";

    cout << "\n\n";
    globalK.SortIt();
    //SortCOO(globalK.get_x(), globalK.get_y(), globalK.get_data(), loads.get_size(), globalK.get_size());


    end_time = clock();
    cout<< "Time(SortCOO): "<< end_time - start_time<< " ms\n";

    ApplyConstraints(globalK, constraints, loads.get_size());
    cout << "ApplyConstraints success\n";

    end_time = clock();
    cout<< "Time(ApplyConstraints): "<< end_time - start_time<< " ms\n";


    int nonzero = globalK.CountNonZero();
    cout << "nonzero = " << nonzero << endl;
    SparseMatrixCOO globalK2(nonzero);
    globalK2 = globalK.DeleteZeros();


    globalK2.CGM_solve(loads, displacements, loads.get_size(), 1e-10);

    end_time = clock();
    cout<< "Time(SolveSparseSystem): "<< end_time - start_time<< " ms\n";


    //displacements.Show();


    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<float> sigma_mises;
    std::vector<float> epsilon_mises;

    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, D, elements, displacements);

    MakeVTKfile(output_vtk, nodesX, nodesY, nodesZ, elements, displacements, Stress, sigma_mises, Deformation, epsilon_mises);


    return 0;
}
