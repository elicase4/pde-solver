#ifndef PDESOLVER_SOLVER_CONFIG_LOGGINGCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_LOGGINGCONFIG_HPP

#include <string>

namespace pdesolver {
	namespace solver {
		namespace config {

			struct LoggerConfig {

				enum class Type {
					Console,
					None
				};

				Type type = Type::Console;

			}; // struct LoggerConfig

			struct LoggingConfig {

				LoggerConfig solver;
				LoggerConfig driver;

				// mirror console output to this file too; empty = no file mirroring
				std::string outputFile;

			}; // struct LoggingConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
