#ifndef APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <vector>

#include "application/heateq/config/ConductivityConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class ConductivityConfigParser {
				public:

					static ConductivityConfig::Type parseConductivityType(const std::string& str);

					static ConductivityConfig parse(const YAML::Node& node);
				
				}; // class ConductivityConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
