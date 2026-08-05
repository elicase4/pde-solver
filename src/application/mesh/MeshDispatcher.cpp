#include "application/mesh/MeshDispatcher.hpp"
#include "application/mesh/MeshConfig.hpp"

#include "application/mesh/GmshMeshGenerator.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "io/MeshIO.hpp"

#include <iostream>
#include <stdexcept>

namespace pdesolver {
	namespace application {
		namespace mesh {

			bool MeshDispatcher::run(const MeshConfig& config) {

				pdesolver::mesh::Mesh mesh;

				switch (config.type) {

					case MeshConfig::Type::Block2D: {

						const auto& b = config.block2D;

						pdesolver::mesh::generator::BlockMesh2D gen{b.nx, b.ny, b.xmin, b.xmax, b.ymin, b.ymax, b.Px, b.Py};
						mesh = gen.generate();

						std::cout << "[mesh] Block2D: " << b.nx << "x" << b.ny << " elements generated\n";
						break;
					}

					case MeshConfig::Type::Block3D: {
						// TODO: add BlockMesh3D generator and wire up here.
						throw std::runtime_error("MeshDispatcher: Block3D generator not yet implemented");
					}

					case MeshConfig::Type::Gmsh: {

						if (config.inputFile.empty()) {
							throw std::runtime_error("MeshDispatcher: Gmsh type requires 'file' field in config");
						}

						// TODO: add physical group mapping argument
						pdesolver::application::mesh::GmshMeshGenerator gen{config.inputFile};
						mesh = gen.generate();

						std::cout << "[mesh] Gmsh import: " << config.inputFile << " → " << mesh.data.numElements << " elements\n";
						break;
					}

					default:
						throw std::runtime_error("MeshDispatcher: unknown mesh type");
				}

				if (!mesh.isValid()) {
					std::cerr << "[mesh] error: generated mesh failed validity check\n";
					return false;
				}

				io::MeshIO::writeBinary(mesh, config.outputFile);

				std::cout << "[mesh] wrote " << config.outputFile << "\n";

				return true;
			}

		} // namespace mesh
	} // namespace application
} // namespace pdesolver
