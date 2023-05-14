cmake_minimum_required (VERSION 3.5)

project(segyio2)

set(SEGYIO2_PUBLIC_HEADERS
	include/segy.h
	include/util.h
)

set(SEGYIO2_SOURCES
	src/segy.c
)

add_library(segyio2 ${SEGYIO2_PUBLIC_HEADERS} ${SEGYIO2_SOURCES})
target_include_directories(segyio2 PUBLIC include PRIVATE include src)
