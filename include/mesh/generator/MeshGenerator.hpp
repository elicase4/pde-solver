#ifndef PDESOLVER_MESH_GENERATOR_MESHGENERATOR_HPP
#define PDESOLVER_MESH_GENERATOR_MESHGENERATOR_HPP

#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace mesh {
		namespace generator {

			// Runs once per mesh (not a hot-loop path), so virtual dispatch here is fine --
			// same reasoning already applied to LinearSolverRunner.
			class MeshGenerator {
			public:

				virtual ~MeshGenerator() = default;

				// physicalGroupMap only matters to import-based generators (e.g. Gmsh); it's
				// part of the shared interface so MeshDispatcher can call generate() uniformly
				// across generator types. Parametric generators simply ignore it.
				virtual Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const = 0;

			}; // class MeshGenerator

		} // namespace generator
	} // namespace mesh
} // namespace pdesolver

#endif
