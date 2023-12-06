cmake_minimum_required(VERSION 3.10)

project(fem_utils LANGUAGES C CXX) #CUDA

set(SOURCES
    include/fem_utils_cpu.h
    include/fem_utils.h
    src/fem_utils.cpp
    src/cpu_fem_utils.cpp
)

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set(LINKED_LIBS
    matrix_pack
    datakeeper
)

add_library(fem_utils)

if (WITH_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DWITH_CUDA)

    cuda_compile(cuda_objs include/cu/fem_utils_kernels.cu)
    list(APPEND SOURCES ${cuda_objs})

    list(APPEND SOURCES
        include/cu/fem_utils_kernels.h
        include/cu/cuda_fem_utils.h
        src/cu/cuda_fem_utils.cpp
    )

    list(APPEND INCLUDE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include/cu
            ${CUDA_INCLUDE_DIRS}
    )

    list(APPEND LINKED_LIBS
        ${CUDA_LIBRARIES}
    )

    set_target_properties(fem_utils PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

message(STATUS "Sources: ${SOURCES}")
message(STATUS "CUDA_LIBRARIES: ${CUDA_LIBRARIES}")

target_include_directories(fem_utils PUBLIC ${INCLUDE_DIRS})
target_sources(fem_utils PRIVATE ${SOURCES})
target_link_libraries(fem_utils PUBLIC ${LINKED_LIBS})

if (WITH_TESTS)
    add_executable(test_common_fem_utils test/test_fem_utils.cpp)
    target_link_libraries(test_common_fem_utils fem_utils gtest femfunc)
    add_test(NAME test_common_fem_utils COMMAND test_common_fem_utils)
endif()
