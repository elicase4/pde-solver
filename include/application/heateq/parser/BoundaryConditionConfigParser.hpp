#ifndef PDESOLVER_APPLICATION_HEATEQ_PARSER_BOUNDARYCONDITIONCONFIGPARSER_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PARSER_BOUNDARYCONDITIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include "application/heateq/config/BoundaryConditionConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace parser {

				class BoundaryConditionConfigParser {
				public:

					static BoundaryConditionConfig::Type parseBoundaryConditionType(const std::string& str);
					
					static BoundaryConditionConfig::Type parseBoundaryConditionForm(const std::string& str);
					
					static BoundaryConditionConfig parse(const YAML::Node& node);
					
				}; // class BoundaryConditionConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
