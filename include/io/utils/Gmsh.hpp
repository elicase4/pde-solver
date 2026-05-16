#ifndef PDESOLVER_IO_UTILS_GMSH_HPP
#define PDESOLVER_IO_UTILS_GMSH_HPP

#include <vector>

#include "mesh/exchange/gmsh/ElementType.hpp"

namespace pdesolver {
	namespace io {
		namespace gmsh {

			mesh::exchange::gmsh::ElementType elementTypeFromGmsh(int type);

			Index nodesPerElement(mesh::exchange::gmsh::ElementType type);

			Index parametricDimension(mesh::exchange::gmsh::ElementType type);

			Index facesPerElement(mesh::exchange::gmsh::ElementType type);

			std::vector<Index> basisOrder(mesh::exchange::gmsh::ElementType type);

			std::vector<Index> reorderToSolver(Index* conn, mesh::exchange::gmsh::ElementType type);

		} // namespace gmsh
	} // namespace io
} // namespace pdesolver

#endif
