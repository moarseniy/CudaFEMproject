cmake_minimum_required(VERSION 3.5)

project(FemProject)

option(TOOLS_TIMER ON)
option(FEM3D ON)

if (TOOLS_TIMER)
    add_definitions(-DTOOLS_TIMER)
endif()

if (FEM3D)
    add_definitions(-DDIM=3)
else ()
    add_definitions(-DDIM=2)
endif()

add_subdirectory(Linal2)
add_subdirectory(fem)
add_subdirectory(VTKfile)
add_subdirectory(femfunc)
add_subdirectory(datakeeper)
add_subdirectory(Tools)
add_subdirectory(CudaLib)
add_subdirectory(femstruct)
