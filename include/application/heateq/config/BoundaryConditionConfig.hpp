#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"
#include "application/heateq/config/ConductivityConfig.hpp"
#include "solver/config/NodalFieldReadConfig.hpp"
#include "solver/config/NodalFieldWriteConfig.hpp"

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

					enum class Form {
						FluxBC,
						ValueBC
					}; // enum class Form

					Int boundaryID;

					Type type;

					// reuses solver::config::NodalFieldReadConfig::Mode for consistency with
					// IC/source, but not the struct itself -- Flux needs a per-component vector
					// expression, which doesn't fit NodalFieldReadConfig's scalar shape.
					solver::config::NodalFieldReadConfig::Mode mode = solver::config::NodalFieldReadConfig::Mode::Expression;

					// mode == Expression: Value takes a single scalar expression, Flux takes
					// one expression per spatial component. mode == File: file, regardless of type.
					std::string expression;
					std::vector<std::string> fluxExpression;
					std::string file;

					// optional one-shot export of the resolved BC value -- Value only; Flux is a
					// shape mismatch (SpatialDim-wide, not NumDOFs-wide), see HeatProblem.tpp.
					solver::config::NodalFieldWriteConfig write;

					std::vector<Form> forms;

					ConductivityConfig::Type model;

				}; // struct BoundaryConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
