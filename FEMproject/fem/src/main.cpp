
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>

#include "Tools.h"
#include "Linal2.h"
//#include <cmath>
#include "femfunc.h"
#include "VTKfile.h"
#include "init.h"

using namespace std;

int main(void) {
    #ifdef TOOLS_TIMER
        Timer timer(__func__);
    #endif

    callCudaKernel();

    std::cout << DIM << std::endl;

    std::string name = "bulk_test";

    std::string directory = "C:/Users/mokin/Desktop/git/CudaFEMproject/prepared_meshes/";
    std::string output_vtk = "C:/Users/mokin/Desktop/git/CudaFEMproject/final_results/results.vtk";
    std::string output_results = "C:/Users/mokin/Desktop/git/CudaFEMproject/final_results/output.txt";

    float poissonRatio = 0.3, youngModulus = 2e+11;

    FEMdataKeeper FEMdata;
    FEMdata.ParseFiles(directory, name, poissonRatio, youngModulus);
    FEMdata.ShowInfo();

    MyArray displacements(3 * FEMdata.nodesCount);
    CalculateFiniteElementMethod(FEMdata, displacements);

    //displacements.Show();

    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<float> sigma_mises;
    std::vector<float> epsilon_mises;

    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, FEMdata.D, FEMdata.elements, displacements);

    MakeVTKfile(output_vtk, FEMdata.nodesX, FEMdata.nodesY, FEMdata.nodesZ, FEMdata.elements,
                displacements, Stress, sigma_mises, Deformation, epsilon_mises);

    return 0;
}
