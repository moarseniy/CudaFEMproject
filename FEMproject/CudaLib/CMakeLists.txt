cmake_minimum_required(VERSION 3.5)


CUDA_ADD_LIBRARY(
  cuda_lib
  include/init.h
  src/init.cu
)

TARGET_INCLUDE_DIRECTORIES(cuda_lib PUBLIC include "${CUDA_INCLUDE_DIRS}")
TARGET_LINK_LIBRARIES(cuda_lib femfunc matrix_pack cuda_matrix_pack ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} Tools femstruct datakeeper)
