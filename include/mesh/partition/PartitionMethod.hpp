#ifndef PDESOLVER_MESH_PARTITION_PARTITIONMETHOD_HPP
#define PDESOLVER_MESH_PARTITION_PARTITIONMETHOD_HPP

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace mesh {
		namespace partition {

			class PartitionMethod {
			public:
				void RCB(); // reverse coordinate bisection
			}; // class PartitionMethod

		} // namespace partition
	} // namespace mesh
} // namespace pdesolver

#endif
