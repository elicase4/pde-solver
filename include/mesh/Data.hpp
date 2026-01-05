#ifndef PDESOLVER_MESH_DATA_HPP
#define PDESOLVER_MESH_DATA_HPP

#include "core/Types.hpp"
#include <vector>

namespace pdesolver {
	
	namespace mesh {

		struct Data {
			
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
			std::vector<Int> rng;
			std::vector<Real> C;
		
		}; // struct Data
	
	} // namespace mesh

} // namespace pdesolver

#endif
