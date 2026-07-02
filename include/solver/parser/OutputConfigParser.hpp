#ifndef PDESOLVER_SOLVER_PARSER_OUTPUTCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_OUTPUTCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/OutputConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {
	
			class OutputConfigParser {
			public:

				static OutputConfig parse(const YAML::Node& node);

			}; // class OutputConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
