
#include <iostream>
#include <string>

#include <Tools.h>
#include <femfunc.h>
#include <datakeeper.h>

int main(int argc, char *argv[]) {
  CheckRunTime(__func__)

//  if (argc == 1) {
//    std::cout << "RUN SUCCESS\n";
//    return 0;
//  }

  DEVICE_NAME deviceType = CPU;
  if (deviceType == CUDA)
    CUDA_Matrix::_initCUDA();
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config3D_dyn.json";
  std::string task_name = "test_bulk";//"test_rect_pcg"; small_test bulk_pressure test_bulk
  dataKeeper dk(config_path, task_name);

  if (dk.getTaskType() == "static")
    gpuCalculateFEM_EbE_vec(deviceType, dk, true);
  else
    gpuCalculateFEM_DYN(deviceType, dk, true);

  ResultsDataKeeper RESdata2(false, false, false,
                            dk.get_nodesCount());

  WriteResults(dk, RESdata2);

  return 0;
}
