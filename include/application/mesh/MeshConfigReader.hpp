#ifndef PDESOLVER_APPLICATION_MESH_MESHCONFIGREADER_HPP
#define PDESOLVER_APPLICATION_MESH_MESHCONFIGREADER_HPP

#include <string>

#include "application/mesh/MeshConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			class MeshConfigReader {
			public:

				static MeshConfig read(const std::string& filename);

			}; // class MeshConfigReader

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
