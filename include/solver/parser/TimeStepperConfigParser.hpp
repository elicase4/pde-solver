#ifndef PDESOLVER_SOLVER_PARSER_TIMESTEPPERCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_TIMESTEPPERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/TimeStepperConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {

			class TimeStepperConfigParser {
			public:

				static TimeStepperConfig::Type parseTimeStepperType(const std::string& str);

				static TimeStepperConfig parse(const YAML::Node& node);

			}; // class TimeStepperConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
