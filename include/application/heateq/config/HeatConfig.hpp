#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATCONFIG_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATCONFIG_HPP

#include <string>
#include <vector>

#include "application/heateq/config/BoundaryConditionConfig.hpp"
#include "application/heateq/config/ConductivityConfig.hpp"
#include "application/heateq/config/SourceConfig.hpp"

#include "solver/config/DiscretizationConfig.hpp"
#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/config/OutputConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			struct HeatConfig {

				std::string meshFile;

				solver::config::DiscretizationConfig discretization;

				solver::config::LinearSolverConfig linearSolver;
				
				solver::config::TimeStepperConfig transient;

				solver::config::OutputConfig output;

				ConductivityConfig conductivity;

				SourceConfig source;

				std::vector<BoundaryConditionConfig> boundaryConditions;

			}; // struct HeatConfig

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
