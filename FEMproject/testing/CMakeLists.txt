cmake_minimum_required(VERSION 3.5)

ADD_EXECUTABLE(testing src/testing.cpp)
TARGET_LINK_LIBRARIES(testing Linal2 femfunc VTKfile Tools cuda_lib datakeeper)

