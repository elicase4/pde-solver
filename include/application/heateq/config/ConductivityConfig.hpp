#ifndef PDESOLVER_APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct ConductivityConfig {

					enum class Type {
						Default,
						Constant,
						Anisotropic
					}; // enum class Type

					Type type;

					Real value;

					std::vector<std::vector<Real>> tensor;

					std::string unit; // required, SI unit label -- e.g. "W/(m*K)"

				}; // struct ConductivityConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
