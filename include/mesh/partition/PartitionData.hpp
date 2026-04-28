#ifndef PDESOLVER_MESH_PARTITION_PARTITIONDATA_HPP
#define PDESOLVER_MESH_PARTITION_PARTITIONDATA_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace mesh {
		namespace partition {

			struct PartitionData {

				Index numParts; // number of paritions
				Index localPart; // rank parition ID
	
				std::vector<Index> elementOwner; // element -> partition
				std::vector<Index> nodeOwner; // node -> partition

				std::vector<Index> ghostNodes; // nodes owned by neighbor ranks needed on local rank
				std::vector<Index> sharedNodes; // nodes needed on local rank shared with neighbor ranks

			}; // struct PartitionData

		} // namespace partition
	} // namespace mesh
} // namespace pdesolver

#endif
