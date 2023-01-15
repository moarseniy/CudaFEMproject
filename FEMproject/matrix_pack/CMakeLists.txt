project(matrix_pack)

add_library(matrix_pack)
target_include_directories(matrix_pack
    PUBLIC include
    PRIVATE src
)

target_sources(matrix_pack PRIVATE
    include/matrix_pack/matrix_pack.h
    src/matrix_pack.cpp
    src/cpu_matrix.cpp
)

CUDA_ADD_LIBRARY(
    cuda_matrix_pack
    include/cuda_matrix_pack/cuda_matrix.h
    include/cuda_matrix_pack/matrix_kernels.h
    src/cu/matrix_kernels.cu
    src/cu/cuda_matrix.cpp
)

target_include_directories(cuda_matrix_pack
    PUBLIC include/cuda_matrix_pack "${CUDA_INCLUDE_DIRS}"
    PRIVATE src/cu)

target_link_libraries(cuda_matrix_pack matrix_pack ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})

add_executable(test_common_matrix_pack test/test_matrix_pack.cpp)
target_link_libraries(test_common_matrix_pack matrix_pack cuda_matrix_pack gtest)
add_test(NAME test_common_matrix_pack COMMAND test_common_matrix_pack)
