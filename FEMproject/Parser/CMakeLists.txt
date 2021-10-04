cmake_minimum_required (VERSION 3.5)


add_library(
	Parser 
	src/Parser.h 
	src/Parser.cpp
	)

target_include_directories(Parser PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(Parser Linal2 femfunc)
