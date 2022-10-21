cmake_minimum_required (VERSION 3.5)


add_library(
  Tools
  src/Tools.h
  src/Tools.cpp
)

target_include_directories(Tools PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(Tools)
