#include "application/heateq/HeatApplication.hpp"
#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"

#include "solver/logging/LoggerFactory.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	if (argc < 2) {
		std::cerr << "Usage: heateq <config.yaml>\n";
		return EXIT_FAILURE;
	}

	// no LoggingConfig exists until the config actually parses -- plain cerr here is the
	// honest boundary, not an oversight.
	pdesolver::application::heateq::config::HeatConfig cfg;

	try {
		cfg = pdesolver::application::heateq::parser::HeatConfigParser::read(argv[1]);
	} catch (const std::exception& e) {
		std::cerr << "[heateq] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}

	const auto logger = pdesolver::solver::logging::makeDriverLogger(cfg.logging.driver, "heateq");

	try {

		pdesolver::application::heateq::HeatApplication app(cfg);
		return app.run();

	} catch (const std::exception& e) {

		logger.error(e.what());
		return EXIT_FAILURE;

	}

}

pdesolver::application::heateq::HeatApplication::HeatApplication(const pdesolver::application::heateq::config::HeatConfig& config) : config_(config) {}

int pdesolver::application::heateq::HeatApplication::run() {

	bool success = false;

	success = pdesolver::application::heateq::HeatDispatcher::run(config_);

	return success ? EXIT_SUCCESS : EXIT_FAILURE;

}
