
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
#include "datakeeper.h"

using namespace std;

int main(void) {
    CheckRunTime(__func__)

    std::string name = "bulk_test";

    std::string directory = "C:/Users/mokin/Desktop/git/CudaFEMproject/prepared_meshes/";
    std::string output_vtk = "C:/Users/mokin/Desktop/git/CudaFEMproject/final_results/results.vtk";
    std::string output_results = "C:/Users/mokin/Desktop/git/CudaFEMproject/final_results/output.txt";

    float poissonRatio = 0.3, youngModulus = 2e+11;

    FEMdataKeeper FEMdata;
    FEMdata.ParseFiles(directory, name, poissonRatio, youngModulus);
    FEMdata.ShowInfo();

    CalculateFiniteElementMethod(FEMdata);

    MakeResults(FEMdata, output_vtk);

    return 0;
}
