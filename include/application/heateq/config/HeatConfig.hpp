#ifndef PDESOLVER_APPLICATION_HEATEQ_CONFIG_HEATCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_CONFIG_HEATCONFIG_HPP

#include <string>
#include <vector>

#include "application/heateq/config/BoundaryConditionConfig.hpp"
#include "application/heateq/config/ConductivityConfig.hpp"
#include "application/heateq/config/SourceConfig.hpp"

#include "solver/config/MeshConfig.hpp"
#include "solver/config/DiscretizationConfig.hpp"
#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/config/OutputConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace config {

				struct HeatConfig {

					solver::config::mesh mesh;

					solver::config::DiscretizationConfig discretization;

					solver::config::SolverConfig solver;
					
					solver::config::OutputConfig output;

					std::vector<BoundaryConditionConfig> boundaryConditions;

					ConductivityConfig conductivity;
					
					SourceConfig source;

				}; // struct HeatConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
