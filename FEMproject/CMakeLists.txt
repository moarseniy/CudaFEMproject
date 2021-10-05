cmake_minimum_required(VERSION 3.5)

project(FemProject)

if (TOOLS_TIMER)
    add_definitions(-DTOOLS_TIMER)
endif()

add_subdirectory(Linal2)
add_subdirectory(fem)
add_subdirectory(VTKfile)
add_subdirectory(femfunc)
add_subdirectory(Tools)


