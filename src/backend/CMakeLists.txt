add_subdirectory(serial)

if (ENABLE_OMP)
	add_subdirectory(openmp)
endif()

if (ENABLE_OMP)
	add_subdirectory(mpi)
endif()

if (ENABLE_CUDA)
	add_subdirectory(cuda)
endif()
