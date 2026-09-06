#ifndef APPLICATION_HEATEQ_PARSER_INITIALCONDITIONCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_INITIALCONDITIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/InitialConditionConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class InitialConditionConfigParser {
				public:

					static config::InitialConditionConfig parse(const YAML::Node& node);

				}; // class InitialConditionConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
