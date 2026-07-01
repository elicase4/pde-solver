#ifndef APPLICATION_HEATEQ_BOUNDARYCONDITIONCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_BOUNDARYCONDITIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include "application/heateq/BoundaryConditionConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class BoundaryConditionConfigParser {
			public:

				static BoundaryConditionConfig::Type parseBoundaryConditionType(const std::string& str);
				
			}; // class BoundaryConditionConfigParser

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
