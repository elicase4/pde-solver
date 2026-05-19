#ifndef PDESOLVER_IO_UTILS_GMSH_HPP
#define PDESOLVER_IO_UTILS_GMSH_HPP

#include <vector>

#include "mesh/exchange/gmsh/ElementType.hpp"

namespace pdesolver {
	namespace io {
		namespace gmsh {

			struct VecHash {
				
				std::size_t operator()(const std::vector<Index>& v) const noexcept {
					
					std::size_t seed = v.size();
					for (auto x : v){
						seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
					}
					
					return seed;
				}

			}; // struct VecHash

			mesh::exchange::gmsh::ElementType elementTypeFromGmsh(int type);

			Index nodesPerElement(mesh::exchange::gmsh::ElementType type);

			Index parametricDimension(mesh::exchange::gmsh::ElementType type);

			Index facesPerElement(mesh::exchange::gmsh::ElementType type);

			std::vector<Index> basisOrder(mesh::exchange::gmsh::ElementType type);

			std::vector<Index> localFaceNodes(const Index* elemNodes, mesh::exchange::gmsh::ElementType type, Index face);

		} // namespace gmsh
	} // namespace io
} // namespace pdesolver

#endif
