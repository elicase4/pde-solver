#ifndef PDESOLVER_SOLVER_PARSER_NONLINEARSOLVERCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_NONLINEARSOLVERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/NonlinearSolverConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {

			class NonlinearSolverConfigParser {
			public:

				static NonlinearSolverConfig::Type parseNonlinearSolverType(const std::string& str);

				static NonlinearSolverConfig parse(const YAML::Node& node);

			}; // class NonlinearSolverConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
