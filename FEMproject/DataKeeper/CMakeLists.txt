cmake_minimum_required (VERSION 3.5)

add_library(
	datakeeper
	src/datakeeper.h 
	src/datakeeper.cpp
	)

target_include_directories(datakeeper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
TARGET_LINK_LIBRARIES(datakeeper Linal2 Tools femstruct)
