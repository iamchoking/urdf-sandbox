cmake_minimum_required(VERSION 3.10)
project(CS672_2024F_DYN)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

create_executable(pogo    src/pogo_sandbox.cpp)

