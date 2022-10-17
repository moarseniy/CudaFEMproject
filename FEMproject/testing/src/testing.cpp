
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
  if (argc == 1) {
    std::cout << "TESTING RUN SUCCESS\n";
    return 0;
  }

//  std::string name = argv[1];
//  std::string project_directory = argv[2];
//  std::string prepared_meshes_directory = argv[3];
//  std::string results_directory = argv[4];

//  std::string output_vtk = results_directory + "/results.vtk";

//  float poissonRatio = std::stof(argv[5]), youngModulus = std::stof(argv[6]);

//  bool withSmooth = false;
//  bool withMises = false;

//  FEMdataKeeper FEMdata(name, project_directory, prepared_meshes_directory, results_directory);
//  FEMdata.ParseFiles(poissonRatio, youngModulus);

//  bool PRINT_DEBUG_INFO = false;

//  //    CalculateFEM(FEMdata);
//  CalculateFEM_EbE_vec(FEMdata, PRINT_DEBUG_INFO);

//  ResultsDataKeeper RESdata(withSmooth, withMises, FEMdata.nodesCount);

//  MakeResults(FEMdata, RESdata);
//  WriteResults(FEMdata, RESdata, output_vtk, PRINT_DEBUG_INFO);

  return 0;
}
