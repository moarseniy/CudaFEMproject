
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
#include "init.h"

int main(int argc, char *argv[]) {
  CheckRunTime(__func__)

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
    //float dt = 4e-06f, endtime = 0.3f; //0.01f; //0.7f;
    float dt = 0.002f, endtime = 0.7f;
          //dt = 0.002f; //0.00243826f;    // from Fidesys
    float rho = 2400.0f;                // Density. Consider making std::stof(argv[7])
    float damping_alpha = 1e-10f, damping_beta = 1e-10f;    // C = alpha * M + beta * K;
    // beta2 = 0.0 -- explicit scheme (assuming both M and C are diagonal -- make sure to lump mass matrix!)
    // implicit scheme: beta1 >= beta2 >= 1/2
//    float beta1 = 0.55f, beta2 = 0.5f;  // implicit
    float beta1 = 0.5f, beta2 = 0.0f;  // explicit

  bool withSmooth = SMOOTH;
  bool withMises = MISES;

  FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
  FEMdata.ParseFiles(poissonRatio, youngModulus);
  FEMdata.ShowInfo();

  bool PRINT_DEBUG_INFO = true;

  //    CalculateFEM(FEMdata);
  //    CalculateFEM_EbE(FEMdata);
  //    CalculateFEM_EbE_vec(FEMdata);

//  CalculateFEM_EbE_vec_GPU(FEMdata, PRINT_DEBUG_INFO);

//      float endtime = 0.2f;
//      float dx = 0.7810f; // 9task_3         // minimal linear size of an element
//      float Vp = std::sqrtf( ( (youngModulus*(1-poissonRatio)) / ((1+poissonRatio)*(1-2*poissonRatio)) ) / rho ); // P-wave velocity
//      float dt_coef = 0.8f;
//      float dt = dt_coef * std::sqrtf(2.0f)/2.0f * dx / Vp;
//      CalculateFEM_dyn_vec(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
//      gpuCalculateFEM_dyn_relaxation(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1);
//      gpuCalculateFEM_dyn_explicit(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1);

  gpuCalculateFEM_explicit(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, PRINT_DEBUG_INFO);

  ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

  MakeResults(FEMdata, RESdata);
  WriteResults(FEMdata, RESdata, output_vtk, PRINT_DEBUG_INFO);

  return 0;
}
