#include "application/mesh/GmshMeshGenerator.hpp"

#include "io/GmshReader.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"
#include "mesh/exchange/gmsh/MeshConverter.hpp"

pdesolver::mesh::Mesh pdesolver::application::mesh::GmshMeshGenerator::generate(const std::unordered_map<Int, Int>& physicalGroupMap) const {

	pdesolver::mesh::exchange::gmsh::IntermediateMesh intermediateMesh;
	pdesolver::io::GmshReader::read(intermediateMesh, inputFile_);

	pdesolver::mesh::Mesh mesh;
	pdesolver::mesh::exchange::gmsh::MeshConverter::toSolverMesh(mesh, intermediateMesh, physicalGroupMap);

	return mesh;

}
