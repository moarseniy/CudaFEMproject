cmake_minimum_required (VERSION 3.5)

if (FEM3D)
    add_library(
            femfunc
            src/femfunc.h
            src/femfunc.cpp
            )
    target_include_directories(femfunc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
else ()
    add_library(
            femfunc
            src2D/femfunc.h
            src2D/femfunc.cpp
            )
    target_include_directories(femfunc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src2D")
endif()


TARGET_LINK_LIBRARIES(femfunc Linal2 Tools datakeeper femstruct VTKfile)
