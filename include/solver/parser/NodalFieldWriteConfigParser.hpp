#ifndef PDESOLVER_SOLVER_PARSER_NODALFIELDWRITECONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_NODALFIELDWRITECONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include "solver/config/NodalFieldWriteConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {

			class NodalFieldWriteConfigParser {
			public:

				// node is the OWNING section (e.g. 'initial_condition:' or a boundary_conditions
				// entry) -- 'write:' itself is optional; absent means write disabled.
				static config::NodalFieldWriteConfig parse(const YAML::Node& node);

			}; // class NodalFieldWriteConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
