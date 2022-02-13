
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include <direct.h>     // for _mkdir (Visual Studio specific)

#include "Tools.h"
#include "Linal2.h"
//#include <cmath>
#include "femfunc.h"
#include "VTKfile.h"
#include "init.h"
#include "datakeeper.h"
#include "tests.h"

int main(void) {
//    std::vector<Element> elems;
//    std::vector<float*> m_e;
//    std::vector<float*> b_e;
//    float soln[18];
//    test_EbePCG_diag(elems, m_e, b_e, soln);

    CheckRunTime(__func__)

    std::string name = "9task40";

    std::string project_directory = "C:/Users/mokin/Desktop/git/CudaFEMproject/";
    std::string results_directory = project_directory + "final_results/" + std::to_string(DIM) + "D/" + name + "/";
    _mkdir(results_directory.c_str());
    std::string output_vtk = results_directory + "results.vtk";

    float poissonRatio = 0.25, youngModulus = 2e+7;

    FEMdataKeeper FEMdata(name, project_directory);
    FEMdata.ParseFiles(poissonRatio, youngModulus);
    FEMdata.ShowInfo();

    CalculateFiniteElementMethod(FEMdata);

    bool withSmooth = true;
    bool withMises = true;
    ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

    MakeResults(FEMdata, RESdata);
    WriteResults(FEMdata, RESdata, output_vtk);
    return 0;
}
