#include "application/mesh/MeshApplication.hpp"
#include "application/mesh/MeshConfigParser.hpp"
#include "application/mesh/MeshDispatcher.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	if (argc < 2) {
		std::cerr << "Usage: mesh <config.yaml>\n";
		return EXIT_FAILURE;
	}

	try {
		
		const auto cfg = pdesolver::application::mesh::MeshConfigParser::read(argv[1]);
		pdesolver::application::mesh::MeshApplication app(cfg);
		return app.run();
	
	} catch (const std::exception& e) {
		
		std::cerr << "[mesh] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;
	
	}

}

pdesolver::application::mesh::MeshApplication::MeshApplication(const MeshConfig& config) : config_(config) {}

int pdesolver::application::mesh::MeshApplication::run() {
	const bool ok = pdesolver::application::mesh::MeshDispatcher::run(config_);
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
