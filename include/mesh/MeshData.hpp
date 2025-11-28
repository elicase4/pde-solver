#ifndef PDESOLVER_MESHDATA_HPP
#define PDESOLVER_MESHDATA_HPP

#include "core/Types.hpp"
#include <vector>

namespace pdesolver {
	
	namespace mesh {

		struct MeshData {
			
			// dimension
			Index parametricDim;
			Index spatialDim;
			
			// order
			std::vector<Index> basisOrder;
			
			// count
			Index numNodes;
			Index numElements;
			Index facesPerElement;
			Index nodesPerElement;

			// mesh data
			std::vector<Real> xyz;
			std::vector<Index> ien;
			std::vector<Index> rng;
			std::vector<Real> C;
		
		}; // struct MeshData
	
	} // namespace mesh

} // namespace pdesolver

#endif
