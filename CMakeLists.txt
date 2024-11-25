cmake_minimum_required(VERSION 3.10)
project(CS672_2024F_DYN)

# [IMPORTANT] LINK raisim so that other programs (like vscode) can see it!
set(CMAKE_PREFIX_PATH ../raisimLib/raisim/linux)
find_package(raisim CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

create_executable(pogo    src/pogo_sandbox.cpp)
create_executable(pogo1D  src/pogo_1D.cpp)

