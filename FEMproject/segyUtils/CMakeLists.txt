cmake_minimum_required (VERSION 3.5)


add_library(
  segyUtils
  src/segyUtils.h
  src/segyUtils.cpp
)

target_include_directories(segyUtils PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(segyUtils matrix_pack cuda_matrix_pack Tools segyio2)
