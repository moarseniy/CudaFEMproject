cmake_minimum_required (VERSION 3.5)


add_library(
  femstruct
  src/femstruct.h
  src/femstruct.cpp
)
target_include_directories(femstruct PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")


TARGET_LINK_LIBRARIES(femstruct matrix_pack Tools)
