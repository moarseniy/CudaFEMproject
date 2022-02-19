
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include <unordered_map>

#include "Tools.h"
#include "Linal2.h"
#include "femfunc.h"
#include "VTKfile.h"
#include "init.h"
#include "datakeeper.h"
#include "tests.h"

int main(int argc, char *argv[]) {
//    std::vector<Element> elems;
//    std::vector<float*> m_e;
//    std::vector<float*> b_e;
//    float soln[18];
//    test_EbePCG_diag(elems, m_e, b_e, soln);

    CheckRunTime(__func__)

    std::string name = argv[1];
    std::string project_directory = argv[2];
    std::string prepared_meshes_directory = argv[3];
    std::string results_directory = argv[4] + std::to_string(DIM) + "D/" + name + "/";

    std::string output_vtk = results_directory + "results.vtk";

    float poissonRatio = std::stof(argv[5]), youngModulus = std::stof(argv[6]);

    bool withSmooth = SMOOTH;
    bool withMises = MISES;

    FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
    FEMdata.ParseFiles(poissonRatio, youngModulus);
    FEMdata.ShowInfo();

    CalculateFiniteElementMethod(FEMdata);

    ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

    MakeResults(FEMdata, RESdata);
    WriteResults(FEMdata, RESdata, output_vtk);

    return 0;
}
