cmake_minimum_required (VERSION 3.5)

if (FEM3D)
    add_library(
            VTKfile
            src/VTKfile.h
            src/VTKfile.cpp
            )
    target_include_directories(VTKfile PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
else ()
    add_library(
            VTKfile
            src2D/VTKfile.h
            src2D/VTKfile.cpp
            )
    target_include_directories(VTKfile PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src2D")
endif()

target_link_libraries(VTKfile Linal2 femfunc Tools)
