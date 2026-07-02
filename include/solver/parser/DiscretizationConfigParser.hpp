#ifndef PDESOLVER_SOLVER_PARSER_DISCRETIZATIONCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_DISCRETIZATIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/DiscretizationConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {
	
			class DiscretizationConfigParser {
			public:

				static DiscretizationConfig parse(const YAML::Node& node);

			}; // class DiscretizationConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
