#ifndef PDESOLVER_PARALLEL_HALOEXCHANGE_HPP
#define PDESOLVER_PARALLEL_HALOEXCHANGE_HPP

#include <mpi.h>

#include "mesh/partition/PartitionData.hpp"

namespace pdesolver {
	namespace parallel {

		template<typename VectorType>
		class HaloExchange {
		public:
			HaloExchange(const PartitionData& partData, MPI_Comm comm);
			
			// reduce to host
			void reverseExchange(VectorType& vec) const;

			// broadcast ownder to ghosts
			void forwardExchange(VectorType& vec) const;

		private:
			struct SendRecvBuffer {};
			std::vector<SendRecvBuffer> sendBufs, recvBufs;
			MPI_Comm comm_;
		
		}; // class HaloExchange

	} // namespace parallel
} // namespace pdesolver

#endif
