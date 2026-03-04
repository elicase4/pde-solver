cmake_minimum_required(VERSION 3.18)

project(PDESolver
	VERSION 0.1.0
	DESCRIPTION "Modular Finite Element PDE Solver with GPU Acceleration"
	LANGUAGES CXX
)

# Cpp Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# build type and basic compiler summary
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# Build Options
option(BUILD_TESTS "Build unit tests" ON)

# Default compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
	add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Global includes
include_directories(${CMAKE_SOURCE_DIR}/include)

# subdirectories
add_subdirectory(src/mesh)
add_subdirectory(src/io)

# build commands
if(BUILD_TESTS)
	enable_testing()
	find_package(GTest REQUIRED)
	include(GoogleTest)
	add_subdirectory(tests)
endif()

# Summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "  PDE Solver Configuration Summary")
message(STATUS "========================================")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Options:")
message(STATUS "    Tests:      ${BUILD_TESTS}")
message(STATUS "========================================")
message(STATUS "")
