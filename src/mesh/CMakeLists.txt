add_library(pde_mesh
	generator/BlockMesh2D.cpp
)

target_include_directories(pde_mesh
	PUBLIC
		${CMAKE_SOURCE_DIR}/include
)
