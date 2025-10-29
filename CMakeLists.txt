cmake_minimum_required(VERSION 3.10)
project(URDF_SANDBOX)

# [IMPORTANT] LINK raisim so that other programs (like vscode) can see it!
set(CMAKE_PREFIX_PATH ../raisimLib/raisim/linux)
find_package(raisim CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

create_executable(pogo    src/pogo_sandbox.cpp)
create_executable(pogo1D  src/pogo_1D.cpp)
create_executable(pogoSC  src/pogo_showcase.cpp)

create_executable(head    src/head_test.cpp)

create_executable(raipal_ragdoll   src/raipal_ragdoll.cpp)
create_executable(raipal   src/raipal.cpp)

create_executable(minitaur src/minitaur.cpp)
create_executable(minitaur_ragdoll src/minitaur_ragdoll.cpp)