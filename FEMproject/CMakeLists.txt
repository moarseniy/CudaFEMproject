cmake_minimum_required(VERSION 3.5)

project(FemProject)

option(TOOLS_TIMER ON)
option(FEM3D ON)

if (TOOLS_TIMER)
    add_definitions(-DTOOLS_TIMER)
endif()

add_subdirectory(Linal2)
add_subdirectory(fem)
add_subdirectory(VTKfile)
add_subdirectory(femfunc)
#add_subdirectory(femstruct)
add_subdirectory(Tools)

