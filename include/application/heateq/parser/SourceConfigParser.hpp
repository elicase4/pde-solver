#ifndef APPLICATION_HEATEQ_PARSER_SOURCECONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_SOURCECONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/SourceConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class SourceConfigParser {
				public:

					static config::SourceConfig::Type parseSourceType(const std::string& str);

					static config::SourceConfig parse(const YAML::Node& node);
				
				}; // class SourceConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
