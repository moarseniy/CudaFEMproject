cmake_minimum_required(VERSION 3.5)

project(FemProject)

# is required for windows
include_directories(thirdparty)

enable_testing()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake targets")

set(SOURCE_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(BINARY_ROOT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(THIRDPARTY_DIR_NAME thirdparty)
set(THIRDPARTY_SOURCE_DIR ${SOURCE_ROOT_DIR}/${THIRDPARTY_DIR_NAME})
set(THIRDPARTY_BINARY_DIR ${BINARY_ROOT_DIR}/${THIRDPARTY_DIR_NAME})

add_subdirectory(
        ${THIRDPARTY_SOURCE_DIR}/gtest
        ${THIRDPARTY_BINARY_DIR}/gtest
        ${THIRDPARTY_ADDING_EXTRA_ARG}
  )


option(WITH_CUDA "Use CUDA for calculations" ON)
option(WITH_ONEAPI "Use ONEAPI for calculations" OFF)
if (WITH_CUDA)
    add_definitions(-DWITH_CUDA=true)
    add_definitions(-DWITH_ONEAPI=false)


elseif (WITH_ONEAPI)
    add_definitions(-DWITH_CUDA=false)
    add_definitions(-DWITH_ONEAPI=true)
endif()

find_package(CUDA QUIET REQUIRED)

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_52" )

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")



option(TOOLS_TIMER CACHED ON)
option(SMOOTH CACHED OFF)
option(MISES CACHED OFF)
option(STRESS_ALONG_AXIS OFF)

if (STRESS_ALONG_AXIS)
    add_definitions(-DSTRESS_ALONG_AXIS=true)
else()
    add_definitions(-DSTRESS_ALONG_AXIS=false)
endif()

if (SMOOTH)
    add_definitions(-DSMOOTH=true)
else()
    add_definitions(-DSMOOTH=false)
endif()

if (MISES)
    add_definitions(-DMISES=true)
else()
    add_definitions(-DMISES=false)
endif()

if (TOOLS_TIMER)
    add_definitions(-DTOOLS_TIMER)
endif()

#add_subdirectory(Linal2)
add_subdirectory(fem)
add_subdirectory(VTKfile)
add_subdirectory(femfunc)
add_subdirectory(DataKeeper)
add_subdirectory(Tools)
add_subdirectory(CudaLib)
add_subdirectory(femstruct)
#add_subdirectory(testing)
add_subdirectory(matrix_pack)
add_subdirectory(fem_utils)

add_subdirectory(thirdparty/json)
add_subdirectory(thirdparty/segyio)

add_subdirectory(thirdparty/segyio2)
add_subdirectory(segyUtils)
