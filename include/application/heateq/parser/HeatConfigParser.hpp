#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATCONFIGPARSER_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATCONFIGPARSER_HPP

#include <string>

#include "application/heateq/config/HeatConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			class HeatConfigParser {
			public:

				static HeatConfig read(const std::string& filename);

			}; // HeatConfigParser

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
