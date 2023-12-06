cmake_minimum_required(VERSION 3.10)
project(matrix_pack LANGUAGES C CXX) #CUDA

option(WITH_CUDA CACHED ON)

# Set the paths to the source files
set(SOURCES
    include/matrix_pack_cpu.h
    include/matrix_pack.h
    src/matrix_pack.cpp
    src/cpu_matrix.cpp
)

set(INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set(LINKED_LIBS
    Tools
)

# Add the library
add_library(matrix_pack)

if (WITH_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DWITH_CUDA)

    cuda_compile(cuda_objs include/cu/matrix_kernels.cu)
    list(APPEND SOURCES ${cuda_objs})

    list(APPEND SOURCES
        include/cu/matrix_kernels.h
        include/cu/cuda_matrix.h
        src/cu/cuda_matrix.cpp
    )

    list(APPEND INCLUDE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include/cu
            ${CUDA_INCLUDE_DIRS}
    )

    list(APPEND LINKED_LIBS
        ${CUDA_LIBRARIES}
        ${CUDA_CUBLAS_LIBRARIES}
        ${CUDA_cusparse_LIBRARY}
        ${CUDA_cusolver_LIBRARY}
    )

    set_target_properties(matrix_pack PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

message(STATUS "Sources: ${SOURCES}")
message(STATUS "CUDA_LIBRARIES: ${CUDA_LIBRARIES}")

target_include_directories(matrix_pack PUBLIC ${INCLUDE_DIRS})
target_sources(matrix_pack PRIVATE ${SOURCES})
target_link_libraries(matrix_pack PUBLIC ${LINKED_LIBS})

if (WITH_TESTS AND WITH_CUDA)
    add_executable(test_common_matrix_pack test/test_matrix_pack.cpp)
    target_link_libraries(test_common_matrix_pack matrix_pack gtest)
    add_test(NAME test_common_matrix_pack COMMAND test_common_matrix_pack)
endif()
