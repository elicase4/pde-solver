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

				// optional -- plain-text mirror of exactly what's printed to console (or would
				// be, if type is None); empty = no file mirroring
				std::string textFile;

			}; // struct LoggerConfig

			// solver logging additionally supports a structured, per-iteration CSV of residual
			// data -- meaningless for driver-level (discrete-event) logging, hence its own type
			// rather than a field on the shared LoggerConfig.
			struct SolverLoggerConfig : LoggerConfig {

				// optional -- one row per solve iteration (outer_tick, iter, per-DOF abs/rel
				// residuals, flops, elapsed); empty = no CSV
				std::string csvFile;

			}; // struct SolverLoggerConfig

			struct LoggingConfig {

				SolverLoggerConfig solver;
				LoggerConfig driver;

			}; // struct LoggingConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
