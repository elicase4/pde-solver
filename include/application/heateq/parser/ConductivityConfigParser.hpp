#ifndef APPLICATION_HEATEQ_CONDUCTIVITYCONFIGPARSER_HPP
#define APPLICATION_HEATEQ_CONDUCTIVITYCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include "application/heateq/ConductivityConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class ConductivityConfigParser {
			public:

				static ConductivityConfig::Type parseConductivityType(const std::string& str);

			}; // class ConductivityConfigParser

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
