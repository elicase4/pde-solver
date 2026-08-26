#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_INITIALCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_INITIALCONDITIONCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct InitialConditionConfig {

					enum class Type {
						Expression,
						File
					}; // enum class Type

					Type type = Type::Expression;

					std::string expression;

					std::string file;

				}; // struct InitialConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
