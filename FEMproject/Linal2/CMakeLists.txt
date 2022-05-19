cmake_minimum_required (VERSION 3.5)


add_library(
  Linal2
  src/Linal2.h
  src/Linal2.cpp
  )

target_include_directories(Linal2 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
TARGET_LINK_LIBRARIES(Linal2 Tools)
