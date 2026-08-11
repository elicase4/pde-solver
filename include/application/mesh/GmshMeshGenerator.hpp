#ifndef PDESOLVER_APPLICATION_MESH_GMSHMESHGENERATOR_HPP
#define PDESOLVER_APPLICATION_MESH_GMSHMESHGENERATOR_HPP

#include <string>

#include "mesh/generator/MeshGenerator.hpp"

namespace pdesolver {
	namespace application {
		namespace mesh {

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
