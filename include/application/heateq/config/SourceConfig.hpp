#ifndef PDESOLVER_APPLICATION_HEATEQ_PARSER_SOURCECONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PARSER_SOURCECONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct SourceConfig {

					std::string expression;

					enum class Type {
						VolumetricHeatSource
					}; // enum class Type
	
					Type type;

				}; // struct SourceConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
