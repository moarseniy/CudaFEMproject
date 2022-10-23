
#include <iostream>
#include <string>

#include "Tools.h"
#include "femfunc.h"
#include "datakeeper.h"
#include "init.h"


int main(int argc, char *argv[]) {
  CheckRunTime(__func__)

  if (argc == 1) {
    std::cout << "RUN SUCCESS\n";
    return 0;
  }

  // ToDO: Don't send so many parameters into the command prompt
  //       Think how to mage through a config file.
  std::string name = argv[1];
  int DIM = std::stoi(argv[2]);
  std::string project_directory = argv[3];
  std::string prepared_meshes_directory = argv[4];
  std::string results_directory = argv[5];
  std::string output_vtk = results_directory + "/results.vtk";
  float poissonRatio = std::stof(argv[6]), youngModulus = std::stof(argv[7]);
  float rho, damping_alpha, damping_beta, dt, endtime, beta1, beta2;
  bool PRINT_DEBUG_INFO = std::atoi(argv[8]);

  bool withStressAlongAxis = STRESS_ALONG_AXIS;
  bool withSmooth = SMOOTH;
  bool withMises = MISES;

  bool isDYN = ( argc > 9 && std::atoi(argv[9]) );
  if (isDYN) {
    rho = std::stof(argv[10]);
    damping_alpha = std::stof(argv[11]);
    damping_beta = std::stof(argv[12]);
    dt = std::stof(argv[13]);
    endtime = std::stof(argv[14]);
    beta1 = std::stof(argv[15]);
    beta2 = std::stof(argv[16]);
  }

  FEMdataKeeper FEMdata(name, DIM, project_directory, prepared_meshes_directory, results_directory);
  FEMdata.ParseFiles(poissonRatio, youngModulus);
  FEMdata.ShowInfo();

  if (isDYN) {
    gpuCalculateFEM_DYN(FEMdata, rho, damping_alpha, damping_beta, endtime, dt, beta1, beta2, PRINT_DEBUG_INFO);
  } else {
    gpuCalculateFEM_EbE_vec(FEMdata, PRINT_DEBUG_INFO);
  }

  ResultsDataKeeper RESdata(withStressAlongAxis, withSmooth, withMises, FEMdata.nodesCount);

  MakeResults(FEMdata, RESdata);
  WriteResults(FEMdata, RESdata, output_vtk, PRINT_DEBUG_INFO);

  return 0;
}
