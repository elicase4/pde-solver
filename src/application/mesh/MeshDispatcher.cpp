#include "application/mesh/MeshDispatcher.hpp"
#include "application/mesh/MeshConfig.hpp"

#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "mesh/exchange/gmsh/MeshConverter.hpp"
#include "io/GmshReader.hpp"
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

						gen.initializeData();
						gen.generateNodes();
						gen.generateElements();
						gen.generateBoundaryTags();

						mesh = static_cast<pdesolver::mesh::Mesh>(gen);

						std::cout << "[mesh] Block2D: " << b.nx << "x" << b.ny << " element generated\n";
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

						io::GmshReader reader;
						pdesolver::mesh::exchange::gmsh::IntermediateMesh intermediateMesh; 
						reader.read(intermediateMesh, config.inputFile);

						// TODO: add physical group mapping argmuent

						pdesolver::mesh::exchange::gmsh::MeshConverter converter;
						converter.toSolverMesh(mesh, intermediateMesh);

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
