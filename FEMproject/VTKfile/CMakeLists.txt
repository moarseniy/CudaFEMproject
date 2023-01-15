cmake_minimum_required (VERSION 3.5)


add_library(
  VTKfile
  src/VTKfile.h
  src/VTKfile.cpp
)

target_include_directories(VTKfile PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(VTKfile matrix_pack cuda_matrix_pack femfunc Tools)
