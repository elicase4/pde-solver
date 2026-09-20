#ifndef APPLICATION_HEATEQ_PARSER_SPECIFICHEATCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_PARSER_SPECIFICHEATCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/SpecificHeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class SpecificHeatConfigParser {
				public:

					static config::SpecificHeatConfig::Type parseSpecificHeatType(const std::string& str);

					static config::SpecificHeatConfig parse(const YAML::Node& node);

				}; // class SpecificHeatConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
