#ifndef PDESOLVER_SOLVER_PARSER_BACKENDCONFIGPARSER_HPP
#define PDESOLVER_SOLVER_PARSER_BACKENDCONFIGPARSER_HPP

#include <string>

#include <yaml-cpp/yaml.h>

#include "solver/config/BackendConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace parser {

			class BackendConfigParser {
			public:

				static config::BackendConfig::Type parseBackendType(const std::string& str);

				static config::BackendConfig parse(const YAML::Node& root);

			}; // class BackendConfigParser

		} // namespace parser
	} // namespace solver
} // namespace pdesolver

#endif
