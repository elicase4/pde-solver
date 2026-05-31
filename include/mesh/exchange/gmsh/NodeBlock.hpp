#ifndef PDESOLVER_MESH_EXCHANGE_GMSH_NODEBLOCK_HPP
#define PDESOLVER_MESH_EXCHANGE_GMSH_NODEBLOCK_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace mesh {
		namespace exchange {
			namespace gmsh {

				struct NodeBlock {

					Int entityDim = -1;
					Int entityTag = -1;

					std::vector<Index> nodeIDs;

				}; // struct NodeBlock

			} // namespace gmsh
		} // namespace exchange
	} // namespace mesh
} // namespace pdesolver

#endif
