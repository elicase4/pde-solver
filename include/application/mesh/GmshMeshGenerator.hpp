#ifndef PDESOLVER_APPLICATION_MESH_GMSHMESHGENERATOR_HPP
#define PDESOLVER_APPLICATION_MESH_GMSHMESHGENERATOR_HPP

#include <string>

#include "mesh/generator/MeshGenerator.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

			// Thin MeshGenerator adapter over GmshReader (pde_io) + MeshConverter (pde_mesh),
			// so MeshDispatcher can treat Gmsh import uniformly with the parametric generators
			// (Block2D, ...). Lives at the application level rather than in pde_mesh because it
			// needs pde_io -- pde_mesh must never depend on pde_io (pde_io already depends on
			// pde_mesh), so this adapter has to sit in a layer that can see both, same as
			// pde_mesh_app already does for MeshDispatcher.
			class GmshMeshGenerator : public pdesolver::mesh::generator::MeshGenerator {
			public:

				explicit GmshMeshGenerator(std::string inputFile): inputFile_(std::move(inputFile)) {};

				pdesolver::mesh::Mesh generate(const std::unordered_map<Int, Int>& physicalGroupMap = {}) const override;

			private:

				std::string inputFile_;

			}; // class GmshMeshGenerator

		} // namespace mesh
	} // namespace application
} // namespace pdesolver

#endif
