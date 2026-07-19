#ifndef PDESOLVER_SOLVER_PARSER_SOLVERCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_SOLVERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/SolverConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {
	
			class SolverConfigParser {
			public:

				static config::SolverConfig parse(const YAML::Node& node);

			}; // class SolverConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
