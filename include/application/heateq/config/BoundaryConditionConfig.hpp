#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct BoundaryConditionConfig {

					enum class Type {
						Value,
						Flux
					}; // enum class Type

					enum class Form {
						FluxBC,
						ValueBC
					}; // enum class Form

					Int boundaryID;

					Type type;

					std::string expression;

					std::vector<Form> forms;
					
					Model model;

				}; // struct BoundaryConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
