#ifndef PDESOLVER_MESH_MESHCONFIGPARSER_HPP
#define PDESOLVER_MESH_MESHCONFIGPARSER_HPP

#include <stdexcept>
#include <string>

#include "mesh/MeshConfig.hpp"

namespace pdesolver {
	namespace mesh {

		class MeshConfigParser {
		public:

			static MeshConfig::Type parseMeshType(const std::string& str);

		}; // class MeshConfigParser

	} // namespace mesh
} // namespace pdesolver

#endif
