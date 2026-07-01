#ifndef PDESOLVER_APPLICATION_MESH_MESHCONFIGPARSER_HPP
#define PDESOLVER_APPLICATION_MESH_MESHCONFIGPARSER_HPP

#include <string>

#include "application/mesh/MeshConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			class MeshConfigParser {
			public:

				static MeshConfig::Type parseMeshType(const std::string& str);

			}; // class MeshConfigReader

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
