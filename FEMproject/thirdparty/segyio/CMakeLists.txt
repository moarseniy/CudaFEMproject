cmake_minimum_required (VERSION 3.5)


add_library(
	segyio
	src/SegyIn.cpp
	src/SegyIn.h
	src/SegyIO.cpp
	src/SegyIO.h
	src/SegyOut.cpp
	src/SegyOut.h
	)

target_include_directories(segyio PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
