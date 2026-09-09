#ifndef PDESOLVER_SOLVER_LOGGING_LOGGERFACTORY_HPP
#define PDESOLVER_SOLVER_LOGGING_LOGGERFACTORY_HPP

#include <string>
#include <utility>

#include "solver/config/LoggingConfig.hpp"

#include "utils/logging/driver/Logger.hpp"
#include "utils/logging/solver/Logger.hpp"

namespace pdesolver {
	namespace solver {
		namespace logging {

			inline utils::logging::driver::Logger makeDriverLogger(const config::LoggerConfig& cfg, std::string tag) {

				const bool consoleEnabled = (cfg.type == config::LoggerConfig::Type::Console);

				if (consoleEnabled || !cfg.textFile.empty()) {
					return utils::logging::driver::Logger(utils::logging::driver::ConsoleLogger{std::move(tag), consoleEnabled, cfg.textFile});
				}

				return utils::logging::driver::Logger(utils::logging::NullLogger{});

			}

		} // namespace logging
	} // namespace solver
} // namespace pdesolver

#endif
