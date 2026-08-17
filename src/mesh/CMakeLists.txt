add_library(pde_mesh
	Mesh.cpp
	exchange/gmsh/IntermediateMesh.cpp
	exchange/gmsh/Utils.cpp
	exchange/gmsh/MeshConverter.cpp
	generator/BlockMesh2D.cpp
	generator/BlockMesh3D.cpp
)

target_include_directories(pde_mesh
	PUBLIC
		${CMAKE_SOURCE_DIR}/include
)
