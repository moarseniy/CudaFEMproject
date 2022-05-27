
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
  float rho, damping_alpha, damping_beta, dt, endtime, beta1, beta2;
  bool PRINT_DEBUG_INFO = std::atoi(argv[7]);
  bool withSmooth = SMOOTH;
  bool withMises = MISES;
  bool isDYN = ( argc > 8 && std::atoi(argv[8]) );
  if (isDYN) {
    rho = std::stof(argv[9]);
    damping_alpha = std::stof(argv[10]);
    damping_beta = std::stof(argv[11]);
    dt = std::stof(argv[12]);
    endtime = std::stof(argv[13]);
    beta1 = std::stof(argv[14]);
    beta2 = std::stof(argv[15]);
  }

  FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
  FEMdata.ParseFiles(poissonRatio, youngModulus);
  FEMdata.ShowInfo();

  if (isDYN) {
    gpuCalculateFEM_DYN(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
  } else {
    CalculateFEM_EbE_vec_GPU(FEMdata, PRINT_DEBUG_INFO);
  }

  ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

  MakeResults(FEMdata, RESdata);
  WriteResults(FEMdata, RESdata, output_vtk, PRINT_DEBUG_INFO);

  return 0;
}
