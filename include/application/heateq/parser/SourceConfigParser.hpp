#ifndef APPLICATION_HEATEQ_PARSER_SOURCECONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_SOURCECONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include "application/heateq/config/ConductivityConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class SourceConfigParser {
				public:

					static SourceConfig::Type parseSourceType(const std::string& str);

					static SourceConfig parse(const YAML::Node& node);
				
				}; // class SourceConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
