cmake_minimum_required (VERSION 3.5)

if (FEM3D)
    add_library(
            datakeeper
            src/datakeeper.h
            src/datakeeper.cpp
            )
    target_include_directories(datakeeper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
else ()
    add_library(
            datakeeper
            src2D/datakeeper.h
            src2D/datakeeper.cpp
            )
    target_include_directories(datakeeper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src2D")
endif()


TARGET_LINK_LIBRARIES(datakeeper Linal2 Tools femstruct)
