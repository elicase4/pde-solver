#ifndef PDESOLVER_APPLICATION_MESH_MESHAPPLICATION_HPP
#define PDESOLVER_APPLICATION_MESH_MESHAPPLICATION_HPP

#include "application/mesh/MeshConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			class MeshApplication {
			public:

				explicit MeshApplication(const MeshConfig& config);
				
				int run();

			private:

				MeshConfig config_;

			}; // class MeshDispatcher

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
