cmake_minimum_required (VERSION 3.5)

if (FEM3D)
  add_library(
    femstruct
    src/femstruct.h
    src/femstruct.cpp
    )
  target_include_directories(femstruct PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
else ()
  add_library(
    femstruct
    src2D/femstruct.h
    src2D/femstruct.cpp
    )
  target_include_directories(femstruct PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src2D")
endif()

TARGET_LINK_LIBRARIES(femstruct Linal2 Tools)
