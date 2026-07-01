#include "application/mesh/MeshConfigParser.hpp"

pdesolver::application::mesh::MeshConfig::Type pdesolver::application::mesh::MeshConfigParser::parseMeshType(const std::string& str) {

	if (str == "block2d") {
		return pdesolver::application::mesh::MeshConfig::Type::Block2D;
	}

	if (str == "block3d") {
		return pdesolver::application::mesh::MeshConfig::Type::Block3D;
	}

	if (str == "gmsh") {
		return pdesolver::application::mesh::MeshConfig::Type::Gmsh;
	}

	throw std::runtime_error("Unknown mesh type: " + str);

}
