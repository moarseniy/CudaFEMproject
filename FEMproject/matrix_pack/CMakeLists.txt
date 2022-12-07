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
    matrix_pack_cuda
    include/matrix_pack_cuda/gpu_matrix.h
    include/matrix_pack_cuda/matrix_kernels.h
    src/cu/matrix_kernels.cu
    src/cu/gpu_matrix.cpp
)

target_include_directories(matrix_pack_cuda
    PUBLIC include/matrix_pack_cuda "${CUDA_INCLUDE_DIRS}"
    PRIVATE src/cu)

#target_sources(matrix_pack_cuda PRIVATE
#    include/matrix_pack_cuda/gpu_matrix.h
#    src/cu/gpu_matrix.cpp
#    )

#TARGET_INCLUDE_DIRECTORIES(matrix_pack_cuda PUBLIC include/matrix_pack_cuda "${CUDA_INCLUDE_DIRS}")
target_link_libraries(matrix_pack_cuda matrix_pack ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})


#not working for VISUAL STUDIO
#set(CUDA_SEPARABLE_COMPILATION ON)
#find_package(CUDA REQUIRED)
#set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc")
##enable_language(CUDA)

#CUDA_SELECT_NVCC_ARCH_FLAGS(GPU_ARCH "Auto")

#list(APPEND CUDA_NVCC_FLAGS ";-O3;${GPU_ARCH};-lineinfo")

#set(CUDA_LIBS ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_curand_LIBRARY})

#include_directories(${CUDA_INCLUDE_DIRS})

#add_library(matrix_pack)
#target_include_directories(matrix_pack
#    PUBLIC include
#    PRIVATE src)

#string(REGEX MATCH "compute_([0-9]*)" _ "${GPU_ARCH}")
#set(ARCH ${CMAKE_MATCH_1})
#message("CUDA Architecture: ${GPU_ARCH}")
#set_target_properties(matrix_pack PROPERTIES CUDA_ARCHITECTURES "${ARCH}")

#target_sources(matrix_pack PRIVATE
#	  include/matrix_pack/matrix_pack.h
#    src/matrix_pack.cpp
#    src/cpu_matrix.cpp
#    )


#message(STATUS "CUDA libs: ${CUDA_LIBS}")
##enable_language(CUDA)

#target_link_libraries(matrix_pack PUBLIC ${CUDA_LIBS})

#target_compile_options(matrix_pack PUBLIC
#    "$<$<BUILD_INTERFACE:$<COMPILE_LANGUAGE:CUDA>>:-Xcompiler=-DDLL_EXPORTS,-D_MWAITXINTRIN_H_INCLUDED>")
#set_target_properties(matrix_pack PROPERTIES POSITION_INDEPENDENT_CODE ON)
#set_target_properties(matrix_pack PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(matrix_pack PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)


#target_sources(matrix_pack PRIVATE
#    # CUDA sources
#    cu/matrix_kernels.cu
#    # C++ sources
#    src/cu/gpu_matrix.cpp
#    # private headers
#    src/cu/matrix_kernels.h
#    )


add_executable(test_common_matrix_pack test/test_matrix_pack.cpp)
target_link_libraries(test_common_matrix_pack matrix_pack matrix_pack_cuda) #Boost::filesystem Boost::boost gtest)
add_test(test_common_matrix_pack test_common_matrix_pack)
