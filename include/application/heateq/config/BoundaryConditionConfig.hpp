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

					// Essential (Dirichlet) vs natural (flux)
					enum class Type {
						Value,
						Flux
					}; // enum class Type

					// Expression vs file
					enum class Mode {
						Expression,
						File
					}; // enum class Mode

					enum class Form {
						FluxBC,
						ValueBC
					}; // enum class Form

					Int boundaryID;

					Type type;
					Mode mode = Mode::Expression;

					// Mode::Expression: Value BCs take a single scalar expression; Flux BCs take one expression per spatial component
					std::string expression;
					std::vector<std::string> fluxExpression;

					// Mode::File: not yet implemented
					std::string file;

					std::vector<Form> forms;
					
					ConductivityConfig::Type model;

				}; // struct BoundaryConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
