cmake_minimum_required(VERSION 3.10)

project(VEQ_S CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(libVEQ)

add_executable(VEQ_S src/main.cpp)

target_link_libraries(VEQ_S libVEQ)
