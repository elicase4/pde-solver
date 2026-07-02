#include "application/heateq/HeatApplication.hpp"
#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	if (argc < 2) {
		std::cerr << "Usage: heateq <config.yaml>\n";
		return EXIT_FAILURE;
	}

	try {
		
		const pdesolver::application::heateq::config::HeatConfig cfg = pdesolver::application::heateq::parser::HeatConfigParser::read(argv[1]);
		pdesolver::application::heateq::HeatApplication app(cfg);
		
		return app.run();
	
	} catch (const std::exception& e) {
		
		std::cerr << "[heateq] fatal error: " << e.what() << "\n";
		
		return EXIT_FAILURE;
	
	}

}

pdesolver::application::heateq::HeatApplication::HeatApplication(const HeatConfig& config) : config_(config) {}

int pdesolver::application::heateq::HeatApplication::run() {

	bool success = false;

	success = HeatDispatcher::run(config_);

	return success ? EXIT_SUCCESS : EXIT_FAILURE;

}
