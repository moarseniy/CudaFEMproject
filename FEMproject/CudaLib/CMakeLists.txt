cmake_minimum_required(VERSION 3.5)
find_package(CUDA QUIET REQUIRED)
#findCUDAToolkit()

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_52" )

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

CUDA_ADD_LIBRARY(
  cuda_lib
  include/init.h
  src/init.cu
)

TARGET_INCLUDE_DIRECTORIES(cuda_lib PUBLIC include)
TARGET_INCLUDE_DIRECTORIES(cuda_lib PUBLIC "${CUDA_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(cuda_lib femfunc Linal2 ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} Tools femstruct datakeeper)
