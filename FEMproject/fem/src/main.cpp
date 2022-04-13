
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
    CheckRunTime(__func__)

    const int N = 7;
    float A[N] = {1, 3, 3, 3, 2, 2, 0}; // input keys
    float B[N] = {9.f, 8.f, 7.f, 6.f, 5.f, 4.f, 3.f}; // input values
    float B2[7];

    //reductionWithMask(B, A, N, B2, 4);


    if (argc == 1) {
        std::cout << "RUN SUCCESS\n";
        return 0;
    }

    std::string name = argv[1];
    std::string project_directory = argv[2];
    std::string prepared_meshes_directory = argv[3];
    std::string results_directory = argv[4];

    std::string output_vtk = results_directory + "/results.vtk";

    float poissonRatio = std::stof(argv[5]), youngModulus = std::stof(argv[6]);
    float rho = 2400.0f;                // Density. Consider making std::stof(argv[7])
    float damping_alpha = 0.0f, damping_beta = 0.0f;    // C = alpha * M + beta * K;
    // beta2 = 0.0 -- explicit scheme (assuming both M and C are diagonal -- make sure to lump mass matrix!)
    // implicit scheme: beta1 >= beta2 >= 1/2
    float beta1 = 0.55f, beta2 = 0.5f;

    bool withSmooth = SMOOTH;
    bool withMises = MISES;

    FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
    FEMdata.ParseFiles(poissonRatio, youngModulus);
    FEMdata.ShowInfo();

//    CalculateFEM(FEMdata);
//    CalculateFEM_EbE(FEMdata);
//    CalculateFEM_EbE_vec(FEMdata);

    float endtime = 0.14f;
    float dx = 0.7810f; // 9task_3         // minimal linear size of an element
    float Vp = std::sqrtf( ( (youngModulus*(1-poissonRatio)) / ((1+poissonRatio)*(1-2*poissonRatio)) ) / rho ); // P-wave velocity
    float dt_coef = 0.8f;
    float dt = dt_coef * std::sqrtf(2.0f)/2.0f * dx / Vp;
    dt = 0.002f; //0.00243826f;    // from Fidesys
    CalculateFEM_dyn_vec(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2);

    ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

    MakeResults(FEMdata, RESdata);
    WriteResults(FEMdata, RESdata, output_vtk);

    return 0;
}
