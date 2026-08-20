#ifndef PDESOLVER_SOLVER_PARSER_LOGGINGCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_LOGGINGCONFIGPARSER_HPP

#include <string>

#include <yaml-cpp/yaml.h>

#include "solver/config/LoggingConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {

			class LoggingConfigParser {
			public:

				static config::LoggerConfig::Type parseLoggerType(const std::string& str);

				// 'logging:' is optional -- an absent node yields the current implicit
				// default behavior (console logging for solver/driver, no file mirror).
				static config::LoggingConfig parse(const YAML::Node& node);

			}; // class LoggingConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
