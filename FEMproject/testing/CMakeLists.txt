cmake_minimum_required(VERSION 3.5)

add_definitions(-DPRINT_DEBUG_INFO=false)

ADD_EXECUTABLE(testing src/testing.cpp)
TARGET_LINK_LIBRARIES(testing Linal2 femfunc VTKfile Tools cuda_lib datakeeper Tests)

