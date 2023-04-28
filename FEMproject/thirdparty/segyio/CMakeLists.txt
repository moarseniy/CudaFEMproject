cmake_minimum_required (VERSION 3.5)

project(segyio)

set(SEGYIO_PUBLIC_HEADERS
	include/SegyIn.h
	include/SegyIO.h
	include/SegyOut.h
)

set(SEGYIO_SOURCES
	src/SegyIn.cpp
	src/SegyIO.cpp
	src/SegyOut.cpp
)

add_library(segyio ${SEGYIO_PUBLIC_HEADERS} ${SEGYIO_SOURCES})
target_include_directories(segyio PUBLIC include PRIVATE include src)