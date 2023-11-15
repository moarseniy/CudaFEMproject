#include <fem_utils/fem_utils.h>
#include <cuda_fem_utils/cuda_fem_utils.h>

#include <gtest/gtest.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include <unordered_set>
#include <map>

#include <femfunc.h>

#define EPS 1e-3

TEST(staticTask, test_rect_pcg) {
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config.json";
  std::string task_name = "test_rect_pcg";
  dataKeeper dk(config_path, task_name);

  gpuCalculateFEM_EbE_vec2(CUDA, dk, false);

  CPU_Matrix cuda_res;
  dk.get_displacements()->copy(cuda_res);
  dk.get_displacements()->setTo(0.f);

  gpuCalculateFEM_EbE_vec2(CPU, dk, false);

  float *cpu_data = dk.get_displacements()->get_data();
  float *gpu_data = cuda_res.get_data();

  for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
    ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: test_rect_pcg failed!";
  }
}

TEST(staticTask, test_bulk) {
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config3D.json";
  std::string task_name = "test_bulk";
  dataKeeper dk(config_path, task_name);

  gpuCalculateFEM_EbE_vec2(CUDA, dk, false);

  CPU_Matrix cuda_res;
  dk.get_displacements()->copy(cuda_res);
  dk.get_displacements()->setTo(0.f);

  gpuCalculateFEM_EbE_vec2(CPU, dk, false);

  float *cpu_data = dk.get_displacements()->get_data();
  float *gpu_data = cuda_res.get_data();

  for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
    ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: test_bulk failed!";
  }
}

TEST(dynamicTask, test_bulk_dyn) {
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config3D_dyn.json";
  std::string task_name = "test_bulk";
  dataKeeper dk(config_path, task_name);

  gpuCalculateFEM_DYN2(CUDA, dk, false);

  CPU_Matrix cuda_res;
  dk.get_displacements()->copy(cuda_res);
  dk.get_displacements()->setTo(0.f);

  gpuCalculateFEM_DYN2(CPU, dk, false);

  ResultsDataKeeper RESdata2(false, false, false,
                            dk.get_nodesCount());

  WriteResults(dk, RESdata2);

  float *cpu_data = dk.get_displacements()->get_data();
  float *gpu_data = cuda_res.get_data();

  for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
    ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: test_bulk_dyn failed!";
  }
}

//TEST(dynamicTask, quasi1Dwave_taller) {
//  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config_dyn_test.json";
//  std::string task_name = "quasi1Dwave_taller";
//  dataKeeper dk(config_path, task_name);

//  gpuCalculateFEM_DYN2(CUDA, dk, false);

//  CPU_Matrix cuda_res;
//  dk.get_displacements()->copy(cuda_res);
//  dk.get_displacements()->setTo(0.f);

//  gpuCalculateFEM_DYN2(CPU, dk, false);

//  float *cpu_data = dk.get_displacements()->get_data();
//  float *gpu_data = cuda_res.get_data();

//  for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
//    ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: test_rect_pcg failed!";
//  }
//}

// TODO: think how to check Klocals
TEST(UtilsFuncs, CalculateKlocals) {
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config3D.json";
  std::string task_name = "test_bulk";
  dataKeeper dk(config_path, task_name);

  ElementsData *elemsData_cuda = ElementsData::setElementsData(CUDA, dk);
  elemsData_cuda->calculateKlocals();
//  elemsData_cuda->calculateDiag(*elemsData_cuda->get_diagK(), 0.f, 1.f);
  elemsData_cuda->getDiagonalElements(*elemsData_cuda->get_Klocals(), *elemsData_cuda->get_diagK());

  ElementsData *elemsData_cpu = ElementsData::setElementsData(CPU, dk);
  elemsData_cpu->calculateKlocals();
//  elemsData_cpu->calculateDiag(*elemsData_cpu->get_diagK(), 0.f, 1.f);
  elemsData_cpu->getDiagonalElements(*elemsData_cpu->get_Klocals(), *elemsData_cpu->get_diagK());

  {
    // calculateAreas
    CPU_Matrix cuda_res;
    elemsData_cuda->get_elementsAreas()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_elementsAreas()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: calculateAreas failed!";
    }
  }
  {
    // genCoordinates test
    CPU_Matrix cuda_res;
    elemsData_cuda->get_coordinates()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_coordinates()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: genCoordinates failed!";
    }
  }
  {
    // genGradientMatrix test
    CPU_Matrix cuda_res;
    elemsData_cuda->get_Blocals()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_Blocals()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: genGradientMatrix failed!";
    }
  }
//  {
//    // calculateDiag

//    std::cout << elemsData_cpu->get_diagK()->min() << " " << elemsData_cpu->get_diagK()->max() << "\n";
//    std::cout << elemsData_cuda->get_diagK()->min() << " " << elemsData_cuda->get_diagK()->max() << "\n";

//    CPU_Matrix cuda_res;
//    elemsData_cuda->get_diagK()->copy(cuda_res);
//    float *gpu_data = cuda_res.get_data();
//    float *cpu_data = elemsData_cpu->get_diagK()->get_data();
//    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
//      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: calculateDiag failed!";
//    }
//  }
//  {
//    // calculateKlocals test
//    CPU_Matrix cuda_res;
//    elemsData_cuda->get_Klocals()->copy(cuda_res);
//    float *gpu_data = cuda_res.get_data();
//    float *cpu_data = elemsData_cpu->get_Klocals()->get_data();
//    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
//      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: calculateKlocals failed!";
//    }
//  }

}

TEST(UtilsFuncs, CalculateFlocals) {
  std::string config_path = "C:/Users/mokin/Desktop/git/CudaFEMproject/configs/run_config3D.json";
  std::string task_name = "test_bulk";
  dataKeeper dk(config_path, task_name);

  ElementsData *elemsData_cuda = ElementsData::setElementsData(CUDA, dk);
  elemsData_cuda->initFlocals(0.f, dk.getWaveletParams());

  ElementsData *elemsData_cpu = ElementsData::setElementsData(CPU, dk);
  elemsData_cpu->initFlocals(0.f, dk.getWaveletParams());

  {
    // genFCoordinates test
    CPU_Matrix cuda_res;
    elemsData_cuda->get_fcoordinates()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_fcoordinates()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: genFCoordinates failed!";
    }
  }
  {
    // calculateLength test
    CPU_Matrix cuda_res;
    elemsData_cuda->get_bEdgesLengths()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_bEdgesLengths()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: calculateLength failed!";
    }
  }
  {
    // calculateFlocals test
    CPU_Matrix cuda_res;
    elemsData_cuda->get_Flocals()->copy(cuda_res);
    float *gpu_data = cuda_res.get_data();
    float *cpu_data = elemsData_cpu->get_Flocals()->get_data();
    for (size_t i(0); i < cuda_res.get_numElements(); ++i) {
      ASSERT_NEAR(gpu_data[i], cpu_data[i], EPS) << "ERROR: calculateFlocals failed!";
    }
  }
}

int main(int argc, char *argv[]) {
  std::srand(time(0));
  CUDA_Matrix::_initCUDA();

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
