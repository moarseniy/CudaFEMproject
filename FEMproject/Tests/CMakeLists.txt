cmake_minimum_required (VERSION 3.5)


add_library(
  Tests
  src/tests.h
  src/tests.cpp
  )

target_include_directories(Tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
TARGET_LINK_LIBRARIES(Tests femstruct Tools)
