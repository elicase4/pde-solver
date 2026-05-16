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
		case ET::TriP1:
		case ET::QuadP1:
		case ET::TetP1:
		case ET::HexP1:
			return {1};

		case ET::LineP2:
		case ET::TriP2:
		case ET::QuadP2:
		case ET::TetP2:
		case ET::HexP2:
			return {1};

		default:
			return {};
	
	}

}

std::vector<Index> pdesolver::io::gmsh::reorderToSolver(Index* conn, pdesolver::mesh::exchange::gmsh::ElementType type){

	using ET = mesh::exchange::gmsh::ElementType;

	switch (type) {

		case ET::QuadP1:
			return {conn[0], conn[1], conn[3], conn[2]};

		case ET::QuadP2:
			return {conn[0], conn[4], conn[1], conn[7], conn[8], conn[5], conn[3], conn[6], conn[2]};

		case ET::HexP1:
			return {conn[0], conn[1], conn[3], conn[2],
					conn[4], conn[5], conn[7], conn[8]};

		case ET::HexP2:
			return {conn[0], conn[8], conn[1], conn[11], conn[20], conn[9], conn[3], conn[10], conn[2],
			        conn[12], conn[24], conn[13], conn[22], conn[26], conn[23], conn[15], conn[25], conn[14],
			        conn[4], conn[16], conn[5], conn[19], conn[21], conn[17], conn[7], conn[18], conn[6]};

		case ET::TriP1:
		case ET::TetP1:
		case ET::TriP2:
		case ET::TetP2:
			return conn;

		default:
			return conn;
	
	}

}
