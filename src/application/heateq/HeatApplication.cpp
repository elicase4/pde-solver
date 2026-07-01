#include "application/heateq/HeatApplication.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"
#include "application/heateq/HeatDispatcher.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {

	if (argc < 2) {
		std::cerr << "Usage: heateq <config.yaml>\n";
		return EXIT_FAILURE;
	}

	try {
		const auto cfg = pdesolver::application::heateq::HeatConfigReader::read(argv[1]);
		pdesolver::application::heateq::HeatApplication app(cfg);
		return app.run();
	}
	catch (const std::exception& e) {
		std::cerr << "[heateq] fatal error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
}

namespace pdesolver {
	namespace application {
		namespace heateq {

			HeatApplication::HeatApplication(const HeatConfig& config) : config_(config) {}

			int HeatApplication::run() {

				bool success = false;

				if (isTransient()) {
					success = HeatDispatcher::runTransient(config_);
				}
				else {
					success = HeatDispatcher::runSteady(config_);
				}

				return success ? EXIT_SUCCESS : EXIT_FAILURE;
			}

			bool HeatApplication::isTransient() const {
				return config_.transient.tf > config_.transient.t0;
			}

		} // namespace heateq
	} // namespace application
} // namespace pdesolver
