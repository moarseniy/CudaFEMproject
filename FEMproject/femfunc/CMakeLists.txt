cmake_minimum_required (VERSION 3.5)


add_library(
  femfunc
  src/femfunc.h
  src/femfunc.cpp
)

target_include_directories(femfunc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

TARGET_LINK_LIBRARIES(femfunc matrix_pack cuda_matrix_pack Tools datakeeper femstruct VTKfile fem_utils cuda_fem_utils segyUtils)
