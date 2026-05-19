#include "mesh/exchange/gmsh/IntermediateMesh.hpp"

std::vector<const pdesolver::mesh::exchange::gmsh::ExchangeBlock*> pdesolver::mesh::exchange::gmsh::IntermediateMesh::cellBlocks() const {

	std::vector<const pdesolver::mesh::exchange::gmsh::ExchangeBlock*> out;

	for (const auto& eb : elementBlocks) {
		if (eb.entityDim == static_cast<Int>(parametricDim)) {
			out.push_back(&eb);
		}
	}

	return out;

}

std::vector<const pdesolver::mesh::exchange::gmsh::ExchangeBlock*> pdesolver::mesh::exchange::gmsh::IntermediateMesh::boundaryBlocks() const {

	std::vector<const pdesolver::mesh::exchange::gmsh::ExchangeBlock*> out;

	if (parametricDim == 0){
		return out;
	}

	for (const auto& eb : elementBlocks) {
		if (eb.entityDim == static_cast<Int>(parametricDim - 1)) {
			out.push_back(&eb);
		}
	}

	return out;

}
