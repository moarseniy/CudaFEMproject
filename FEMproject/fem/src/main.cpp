#include "Linal2.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include "Tools.h"
//#include <cmath>

#include "femfunc.h"
#include "VTKfile.h"

using namespace std;

int main(void) {
    #ifdef TOOLS_TIMER
        Timer timer(__func__);
    #endif

    std::string name = "bulk_test";

    std::string directory = "C:/Users/mokin/Desktop/FiniteElementMethod3D/prepared_meshes/";
    std::string output_vtk = "C:/Users/mokin/Desktop/FiniteElementMethod3D/final_results/results.vtk";
    std::string output_results = "C:/Users/mokin/Desktop/FiniteElementMethod3D/final_results/output.txt";

    FEMdataKeeper FEMdata;
    FEMdata.ParseFiles(directory, name);
    FEMdata.ShowInfo();

    float poissonRatio = 0.3, youngModulus = 2e+11;
    Matrix D(6, 6);

    D(0, 0) = D(1, 1) = D(2, 2) = 1.0;
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = poissonRatio / (1.0 - poissonRatio);
    D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));

    D.scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));

    ///////////////////////////

    std::vector<Triplet> triplets;
    for (std::vector<Element>::iterator it = FEMdata.elements.begin(); it != FEMdata.elements.end(); ++it) {
        it->CalculateStiffnessMatrix(D, triplets, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ);
    }
    cout << "CalculateStiffnessMatrix success\n";
    cout << "Triplets Size = " << triplets.size() << endl;

    MyArray displacements(FEMdata.loads.get_size());
    SparseMatrixCOO globalK(triplets.size());

    globalK.ConvertTripletToSparse(triplets);
    cout << "new size= "<<globalK.get_size()<<"\n";
    cout << "!!!" << sizeof (globalK) << endl;

    globalK.resize();

    globalK.SortIt();
    //SortCOO(globalK.get_x(), globalK.get_y(), globalK.get_data(), loads.get_size(), globalK.get_size());

    ApplyConstraints(globalK, FEMdata.constraints, FEMdata.loads.get_size());


    int nonzero = globalK.CountNonZero();
    cout << "nonzero = " << nonzero << endl;
    SparseMatrixCOO globalK2(nonzero);
    globalK2 = globalK.DeleteZeros();


    globalK2.CGM_solve(FEMdata.loads, displacements, FEMdata.loads.get_size(), 1e-10);

    //displacements.Show();

    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<float> sigma_mises;
    std::vector<float> epsilon_mises;

    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, D, FEMdata.elements, displacements);

    MakeVTKfile(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ, FEMdata.elements,
                displacements, Stress, sigma_mises, Deformation, epsilon_mises);


    return 0;
}
