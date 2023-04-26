project(matrix_pack)

add_library(fem_utils)
target_include_directories(fem_utils
    PUBLIC include
    PRIVATE src
)

target_sources(fem_utils PRIVATE
    include/fem_utils/fem_utils.h
    src/fem_utils.cpp
    src/cpu_fem_utils.cpp
)

target_link_libraries(fem_utils matrix_pack datakeeper)

CUDA_ADD_LIBRARY(
    cuda_fem_utils
    include/cuda_fem_utils/cuda_fem_utils.h
    include/cuda_fem_utils/fem_utils_kernels.h
    src/cu/fem_utils_kernels.cu
    src/cu/cuda_fem_utils.cpp
)

target_include_directories(cuda_fem_utils
    PUBLIC include/cuda_fem_utils "${CUDA_INCLUDE_DIRS}"
    PRIVATE src/cu)

target_link_libraries(cuda_fem_utils fem_utils ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})

add_executable(test_common_fem_utils test/test_fem_utils.cpp)
target_link_libraries(test_common_fem_utils fem_utils cuda_fem_utils gtest femfunc)
add_test(NAME test_common_fem_utils COMMAND test_common_fem_utils)
