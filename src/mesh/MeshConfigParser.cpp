#include "mesh/MeshConfigParser.hpp"

pdesolver::mesh::MeshConfig::Type pdesolver::mesh::MeshConfigParser::parseMeshType(const std::string& str) {

	if (str == "block2d") {
		return pdesolver::mesh::MeshConfig::Type::Block2D;
	}

	if (str == "block3d") {
		return pdesolver::mesh::MeshConfig::Type::Block3D;
	}

	if (str == "gmsh") {
		return pdesolver::mesh::MeshConfig::Type::Gmsh;
	}

	throw std::runtime_error("Unknown mesh type: " + str);

}
