#ifndef PDESOLVER_UTILS_LOGGING_BANNER_HPP
#define PDESOLVER_UTILS_LOGGING_BANNER_HPP

#include <iostream>
#include <string>

namespace pdesolver {
	namespace utils {
		namespace logging {

			// printed once at app launch, before config parsing -- unconditional (no
			// LoggingConfig exists yet to gate it on), same reasoning as the "no LoggingConfig
			// exists until the config actually parses" plain-cerr boundary in main().
			inline void printStartupBanner(const std::string& appName, const std::string& description) {

				const int width = 60;

				std::cout << "\n";
				std::cout << std::string(width, '=') << "\n";
				std::cout << "  Residuum -- " << appName << "\n";
				std::cout << "  " << description << "\n";
				std::cout << std::string(width, '=') << "\n";
				std::cout << "\n";

			}

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
