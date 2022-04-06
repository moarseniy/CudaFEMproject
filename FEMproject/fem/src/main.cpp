
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
    std::string results_directory = argv[4];

    std::string output_vtk = results_directory + "/results.vtk";

    float poissonRatio = std::stof(argv[5]), youngModulus = std::stof(argv[6]);
    float rho = 2400.0f;                // Density. Consider making std::stof(argv[7])
    float damping_alpha = 0.0f, damping_beta = 0.0f;    // C = alpha * M + beta * K;

    bool withSmooth = SMOOTH;
    bool withMises = MISES;

    FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
    FEMdata.ParseFiles(poissonRatio, youngModulus);
    FEMdata.ShowInfo();

//    CalculateFEM(FEMdata);
//    CalculateFEM_EbE(FEMdata);
    CalculateFEM_EbE_vec(FEMdata);

//    float endtime = 0.21f;
//    float dx = 0.7810f; // 9task_3         // minimal linear size of an element
//    float Vp = std::sqrtf( ( (youngModulus*(1-poissonRatio)) / ((1+poissonRatio)*(1-2*poissonRatio)) ) / rho ); // P-wave velocity
//    float dt_coef = 0.8f;
//    float dt = dt_coef * std::sqrtf(2.0f)/2.0f * dx / Vp;
//    dt = 0.01f; //0.00243826f;    // from Fidesys
//    CalculateFEM_dyn(FEMdata, rho, damping_alpha, damping_beta, endtime, dt);

    ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

    MakeResults(FEMdata, RESdata);
    WriteResults(FEMdata, RESdata, output_vtk);

    return 0;
}
