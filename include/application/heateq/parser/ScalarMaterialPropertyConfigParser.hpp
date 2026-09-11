#ifndef APPLICATION_HEATEQ_PARSER_SCALARMATERIALPROPERTYCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_SCALARMATERIALPROPERTYCONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include "application/heateq/config/ScalarMaterialPropertyConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class ScalarMaterialPropertyConfigParser {
				public:

					static config::ScalarMaterialPropertyConfig parse(const YAML::Node& node);

				}; // class ScalarMaterialPropertyConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
