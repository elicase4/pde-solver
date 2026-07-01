#ifndef PDESOLVER_APPLICATION_MESH_MESHDISPATCHER_HPP
#define PDESOLVER_APPLICATION_MESH_MESHDISPATCHER_HPP

#include "application/mesh/MeshConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			class MeshDispatcher {
			public:

				static bool run(const MeshConfig& config);

			}; // class MeshDispatcher

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
