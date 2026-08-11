#ifndef PDESOLVER_MESH_GENERATOR_MESHGENERATOR_HPP
#define PDESOLVER_MESH_GENERATOR_MESHGENERATOR_HPP

#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace mesh {
		namespace generator {

			class MeshGenerator {
			public:

				virtual ~MeshGenerator() = default;

				virtual Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const = 0;

			}; // class MeshGenerator

		} // namespace generator
	} // namespace mesh
} // namespace pdesolver

#endif
