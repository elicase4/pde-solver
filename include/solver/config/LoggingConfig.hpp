#ifndef RESIDUUM_SOLVER_CONFIG_LOGGINGCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_LOGGINGCONFIG_HPP

#include <string>

namespace residuum {
	namespace solver {
		namespace config {

			struct LoggerConfig {

				enum class Type {
					Console,
					None
				};

				Type type = Type::Console;

				// plain-text mirror of the console output; empty means no file mirroring
				std::string textFile;

			}; // struct LoggerConfig

			// adds a structured per-iteration residual CSV, meaningless at driver level
			struct SolverLoggerConfig : LoggerConfig {

				// one row per solve iteration; empty means no CSV
				std::string csvFile;

			}; // struct SolverLoggerConfig

			// one row per accepted step, not per solve attempt
			struct TimeStepperLoggerConfig : LoggerConfig {

				// one row per accepted step; empty means no CSV
				std::string csvFile;

			}; // struct TimeStepperLoggerConfig

			// one row per outer Newton/Picard iteration
			struct NonlinearLoggerConfig : LoggerConfig {

				// one row per outer iteration; empty means no CSV
				std::string csvFile;

			}; // struct NonlinearLoggerConfig

			struct LoggingConfig {

				SolverLoggerConfig solver;
				TimeStepperLoggerConfig timestepper;
				NonlinearLoggerConfig nonlinear;
				LoggerConfig driver;

			}; // struct LoggingConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
