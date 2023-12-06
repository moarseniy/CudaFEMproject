cmake_minimum_required(VERSION 3.5)

ADD_EXECUTABLE(main src/main.cpp)
TARGET_LINK_LIBRARIES(main matrix_pack femfunc Tools datakeeper)

