cmake_minimum_required(VERSION 3.23)
project(simple_formula)

set(CMAKE_CXX_STANDARD 23)
include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
add_executable(simple_formula main.cpp)
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll)
file(COPY ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll DESTINATION ${PROJECT_BINARY_DIR})
