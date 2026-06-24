cmake_minimum_required(VERSION 3.10)
project(URDF_SANDBOX)

set(RAISIM1_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../raisimLib/raisim/linux" CACHE PATH "RaiSim 1.1.8 install prefix")
set(RAISIM2_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../raisim2Lib/raisim" CACHE PATH "RaiSim2 install prefix")
set(URDF_SANDBOX_RAISIM_VERSION "rs1" CACHE STRING "RaiSim backend to build: rs1 or rs2")
set_property(CACHE URDF_SANDBOX_RAISIM_VERSION PROPERTY STRINGS rs1 rs2)

string(TOLOWER "${URDF_SANDBOX_RAISIM_VERSION}" URDF_SANDBOX_RAISIM_VERSION_NORMALIZED)
set(URDF_SANDBOX_RAISIM_TARGET)
set(RAISIM_RAIPAL_USE_RAISIM2 OFF)

# raipal_kinematics is installed with the RaiSim 1 workspace, even when only
# the RaiSim2 source set is being configured.
list(PREPEND CMAKE_PREFIX_PATH "${RAISIM1_PREFIX}")
if(URDF_SANDBOX_RAISIM_VERSION_NORMALIZED MATCHES "^(rs1|raisim1|1)$")
    find_package(raisim 1.1.8 CONFIG REQUIRED)
    set(URDF_SANDBOX_RAISIM_TARGET raisim::raisim)
elseif(URDF_SANDBOX_RAISIM_VERSION_NORMALIZED MATCHES "^(rs2|raisim2|2)$")
    set(URDF_SANDBOX_RAISIM_TARGET raisim2::raisim)
    set(RAISIM_RAIPAL_USE_RAISIM2 ON)
else()
    message(FATAL_ERROR "URDF_SANDBOX_RAISIM_VERSION must be one of: rs1, rs2")
endif()
message(STATUS "URDF sandbox RaiSim backend: ${URDF_SANDBOX_RAISIM_VERSION_NORMALIZED}")
find_package(raipal_kinematics CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(cmake/helper.cmake)

function(add_raisim2_imported_library target library_path)
    if(NOT TARGET "${target}")
        add_library(${target} SHARED IMPORTED GLOBAL)
        set_target_properties(${target} PROPERTIES
            IMPORTED_LOCATION "${library_path}"
            INTERFACE_INCLUDE_DIRECTORIES "${RAISIM2_PREFIX}/include")
    endif()
endfunction()

if(RAISIM_RAIPAL_USE_RAISIM2)
    add_raisim2_imported_library(raisim2::raisimZ "${RAISIM2_PREFIX}/lib/libraisimZ.so")
    add_raisim2_imported_library(raisim2::raisimPng "${RAISIM2_PREFIX}/lib/libraisimPng.so")
    add_raisim2_imported_library(raisim2::raisimMine "${RAISIM2_PREFIX}/lib/libraisimMine.so")
    add_raisim2_imported_library(raisim2::assimp "${RAISIM2_PREFIX}/lib/libassimp.so")
    add_raisim2_imported_library(raisim2::raisim_contact "${RAISIM2_PREFIX}/lib/libraisim_contact.so")
    add_raisim2_imported_library(raisim2::raisim "${RAISIM2_PREFIX}/lib/libraisim.so")

    set_target_properties(raisim2::assimp PROPERTIES
        INTERFACE_LINK_LIBRARIES "${RAISIM2_PREFIX}/lib/libzlibstatic.a;rt")
    set_target_properties(raisim2::raisimPng PROPERTIES
        INTERFACE_LINK_LIBRARIES "raisim2::raisimZ")
    set_target_properties(raisim2::raisim_contact PROPERTIES
        INTERFACE_COMPILE_FEATURES cxx_std_20
        INTERFACE_LINK_LIBRARIES "Eigen3::Eigen")
    set_target_properties(raisim2::raisim PROPERTIES
        INTERFACE_COMPILE_FEATURES cxx_std_20
        INTERFACE_LINK_LIBRARIES "Eigen3::Eigen;raisim2::raisim_contact;raisim2::raisimPng;raisim2::raisimMine;raisim2::assimp")
endif()

set(RAISIM_RAIPAL_TARGET_NAME raisimRaipal)
set(RAISIM_RAIPAL_RAISIM_TARGET ${URDF_SANDBOX_RAISIM_TARGET})
add_subdirectory(
    thirdParty/raisimRaipal
    "${CMAKE_CURRENT_BINARY_DIR}/thirdParty/raisimRaipal"
    EXCLUDE_FROM_ALL)
unset(RAISIM_RAIPAL_TARGET_NAME)
unset(RAISIM_RAIPAL_RAISIM_TARGET)

# Automatically create executables for all .cpp files in src directory and subdirectories
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS "src/*.cpp")

foreach(SOURCE_FILE ${SOURCE_FILES})
    # Get the relative path from src/ and use it as the executable name (replace .cpp with nothing)
    file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src" "${SOURCE_FILE}")
    string(REPLACE ".cpp" "" EXEC_NAME "${REL_PATH}")
    string(REPLACE "/" "_" EXEC_NAME "${EXEC_NAME}")

    create_executable(${EXEC_NAME} ${SOURCE_FILE} ${URDF_SANDBOX_RAISIM_TARGET})
    if(EXEC_NAME MATCHES "^raipal_")
        target_link_libraries(${EXEC_NAME} PRIVATE raipal_kinematics raisimRaipal)
    elseif(EXEC_NAME MATCHES "raipal")
        target_link_libraries(${EXEC_NAME} PRIVATE raipal_kinematics)
    endif()
endforeach()
