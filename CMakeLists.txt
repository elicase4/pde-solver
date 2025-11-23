cmake_minimum_required(VERSION 3.18)

project(FDPDESolver
	VERSION 0.0.1
	DESCRIPTION "Modular Finite Element PDE Solver"
	LANGUAGES CXX
)

# Cpp Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()

# build type and basic compiler summary
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# Build Options
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_BENCHMARKS "Build perforamnce benchmarks" OFF)
option(BUILD_EXAMPLES "Build example applications" ON)
option(ENABLE_OMP "Enable OpenMP backend" OFF)
option(ENABLE_MPI "Enable MPI backend" OFF)
option(ENABLE_CUDA "Enable CUDA backend" OFF)

# Compiler flags
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(compiler_flags OPTIONAL)

# Default compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
	add_compile_options(
		-Wall -Wextra -Wpedantic
		$<$<CONFIG:Debug>:-g -O0>
		$<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
	add_compile_options(
		/W4
		$<$<CONFIG:Debug>:/Od /Zi>
		$<$<CONFIG:Release>:/O2 /DNDEBUG>
	)
endif()

# OpenMP backend
if (ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
	if(OpenMP_CXX_FOUND)
		message(STATUS "OpenMP enabled: ${OpenMP_CXX_VERSION}")
	endif()
endif()

# MPI backend
if (ENABLE_MPI)
	find_package(MPI REQUIRED)
	if(MPI_CXX_FOUND)
		message(STATUS "MPI enabled: ${MPI_CXX_VERSION}")
		message(STATUS "  MPI executable: ${MPIEXEC_EXECUTABLE}")
		message(STATUS "  MPI compile flags: ${MPI_CXX_COMPILE_FLAGS}")
	endif()
endif()

# OpenMP backend
if (ENABLE_CUDA)
	enable_language(CUDA)
	find_package(CUDAToolkit REQUIRED)
	message(STATUS "CUDA enabled: ${CudaToolkit_VERSION}")
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES 75 80 86)
	endif()
	message(STATUS "  CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
endif()

# Test Framework
if (BUILD_TESTS)
	enable_testing()
	
	# Find GoogleTest
	find_package(GTest QUIET)

	if(GTest_FOUND)
		message(STATUS "GoogleTest found: ${GTest_VERSION}")
		include(GoogleTest)
	else()
		message(FATAL_ERROR, "GoogleTest not found and is a requirement.")
	endif()

endif()

# Global includes
include_directories(
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/src
)

# subdirectories
add_subdirectory(src/core)
add_subdirectory(src/mesh)
add_subdirectory(src/fem)
add_subdirectory(src/equations)
add_subdirectory(src/bc)
add_subdirectory(src/linalg)
add_subdirectory(src/utils)
add_subdirectory(src/apps)

# build commands
if(BUILD_TESTS)
	add_subdirectory(tests)
endif()

if(BUILD_BENCHMARKS)
	add_subdirectory(benchmarks)
endif()

if(BUILD_EXAMPLES)
	add_subdirectory(examples)
endif()

# Installation
install(TARGETS
	pde_core pde_mesh pde_fem pde_equations pde_bc pde_linalg pde_utils
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
)

# Headers
install(DIRECTORY include/
	DESTINATION include/pdesolver
	FILES_MATCHING PATTERN "*.hpp" PATTERN ".h"
)

# Summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "  PDE Solver Configuration Summary")
message(STATUS "========================================")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "  Options:")
message(STATUS "    Tests:      ${BUILD_TESTS}")
message(STATUS "    Benchmarks: ${BUILD_BENCHMARKS}")
message(STATUS "    Examples:   ${BUILD_EXAMPLES}")
message(STATUS "    OpenMP:     ${ENABLE_OPENMP}")
message(STATUS "    MPI:        ${ENABLE_MPI}")
message(STATUS "    CUDA:       ${ENABLE_CUDA}")
message(STATUS "========================================")
message(STATUS "")
