cmake_minimum_required(VERSION 3.10)
project(URDF_SANDBOX)

# [IMPORTANT] LINK raisim so that other programs (like vscode) can see it!
set(CMAKE_PREFIX_PATH ../raisimLib/raisim/linux)
find_package(raisim CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
include(cmake/helper.cmake)

# Automatically create executables for all .cpp files in src directory
file(GLOB SOURCE_FILES "src/*.cpp")

foreach(SOURCE_FILE ${SOURCE_FILES})
    # Get filename without path and extension
    get_filename_component(EXEC_NAME ${SOURCE_FILE} NAME_WE)
    
    # Create executable
    create_executable(${EXEC_NAME} ${SOURCE_FILE})
endforeach()