#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"
#include "application/heateq/config/ConductivityConfig.hpp"

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

					// Value BCs: single scalar expression. Flux BCs: one expression per
					// spatial component (BoundaryFluxFunction expects a full vector).
					std::string expression;
					std::vector<std::string> fluxExpression;

					std::vector<Form> forms;
					
					ConductivityConfig::Type model;

				}; // struct BoundaryConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
