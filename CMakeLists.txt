cmake_minimum_required(VERSION 3.23)
project(simple_formula)

set(CMAKE_CXX_STANDARD 23)
add_subdirectory(excerpt)
include_directories(${simple_formula_SOURCE_DIR}/excerpt/inc)
add_executable(simple_formula main.cpp)
target_link_libraries(simple_formula excerpt)
