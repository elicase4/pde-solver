#include "io/utils/Gmsh.hpp"

pdesolver::mesh::exchange::gmsh::ElementType pdesolver::io::gmsh::elementTypeFromGmsh(int type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {
		
		// P1 elements
		case 1:
			return ET::LineP1;
		case 2:
			return ET::TriP1;
		case 3:
			return ET::QuadP1;
		case 4:
			return ET::TetP1;
		case 5:
			return ET::HexP1;
		
		// P2 elements
		case 8:
			return ET::LineP2;
		case 9:
			return ET::TriP2;
		case 10:
			return ET::QuadP2;
		case 11:
			return ET::TetP2;
		case 17:
			return ET::HexP2;

		default:
			return ET::Unknown;
	}

}

Index pdesolver::io::gmsh::nodesPerElement(pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::LineP1:
			return 2;
		case ET::LineP2:
			return 3;

		case ET::TriP1:
			return 3;
		case ET::TriP2:
			return 6;

		case ET::QuadP1:
			return 4;
		case ET::QuadP2:
			return 9;

		case ET::TetP1:
			return 4;
		case ET::TetP2:
			return 10;

		case ET::HexP1:
			return 8;
		case ET::HexP2:
			return 20;

		default:
			return 0;

	}

}

Index pdesolver::io::gmsh::parametricDimension(pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::LineP1:
		case ET::LineP2:
			return 1;

		case ET::TriP1:
		case ET::QuadP1:
		case ET::TriP2:
		case ET::QuadP2:
			return 2;

		case ET::TetP2:
		case ET::HexP2:
		case ET::TetP1:
		case ET::HexP1:
			return 3;

		default:
			return 0;

	}

}

Index pdesolver::io::gmsh::facesPerElement(pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::TriP1:
		case ET::TriP2:
			return 3;

		case ET::QuadP1:
		case ET::QuadP2:
			return 4;

		case ET::TetP1:
		case ET::TetP2:
			return 4;

		case ET::HexP1:
		case ET::HexP2:
			return 6;

		default:
			return 0;

	}

}

std::vector<Index> pdesolver::io::gmsh::basisOrder(pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::LineP1:
			return {1};
		case ET::TriP1:
		case ET::QuadP1:
			return {1,1};
		case ET::TetP1:
		case ET::HexP1:
			return {1,1,1};

		case ET::LineP2:
			return {2};
		case ET::TriP2:
		case ET::QuadP2:
			return {2,2};
		case ET::TetP2:
		case ET::HexP2:
			return {2,2,2};

		default:
			return {};
	
	}

}

std::vector<Index> pdesolver::io::gmsh::localFaceNodes(const Index* elemNodes, pdesolver::mesh::exchange::gmsh::ElementType type, Index face) {

	using ET = pdesolver::mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::TriP1:
			switch (face) {
				case 0: return {elemNodes[0], elemNodes[1]};
				case 1: return {elemNodes[1], elemNodes[2]};
				case 2: return {elemNodes[0], elemNodes[2]};
				default: return {};
			}
		case ET::QuadP1:
			switch (face) {
				case 0: return {elemNodes[0], elemNodes[1]};
				case 1: return {elemNodes[1], elemNodes[3]};
				case 2: return {elemNodes[3], elemNodes[2]};
				case 3: return {elemNodes[2], elemNodes[0]};
				default: return {};
			}
		case ET::TetP1:
			switch (face) {
				case 0: return {elemNodes[0], elemNodes[1], elemNodes[2]};
				case 1: return {elemNodes[0], elemNodes[1], elemNodes[3]};
				case 2: return {elemNodes[1], elemNodes[2], elemNodes[3]};
				case 3: return {elemNodes[0], elemNodes[2], elemNodes[3]};
				default: return {};
			}
		case ET::HexP1:
			switch (face) {
				case 0: return {elemNodes[0], elemNodes[1], elemNodes[2], elemNodes[3]};
				case 1: return {elemNodes[4], elemNodes[5], elemNodes[6], elemNodes[7]};
				case 2: return {elemNodes[0], elemNodes[1], elemNodes[4], elemNodes[5]};
				case 3: return {elemNodes[2], elemNodes[3], elemNodes[6], elemNodes[7]};
				case 4: return {elemNodes[0], elemNodes[2], elemNodes[4], elemNodes[6]};
				case 5: return {elemNodes[1], elemNodes[3], elemNodes[5], elemNodes[7]};
				default: return {};
			}

		case ET::TriP2:
			return {}; // TODO: implement
		case ET::QuadP2:
			return {}; // TODO: implement
		case ET::TetP2:
			return {}; // TODO: implement
		case ET::HexP2:
			return {}; // TODO: implement

		default:
			return {};

}
