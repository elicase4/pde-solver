#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/TimeStepperConfig.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			// Mirrors solver::linear::makeLinearSolverRunner's shape and throw-for-unimplemented
			// pattern. No concrete stepper exists yet -- all four cases throw until #41 lands.
			// TODO: signature will need to grow once a concrete stepper is implemented, to inject
			// whichever of LinearSolverRunner/NonlinearSolverRunner the wrapped stage's resolved
			// SolverMode calls for (see TimeStepperRunner.hpp's composition note).
			template<typename VectorType>
			std::unique_ptr<TimeStepperRunner<VectorType>> makeTimeStepperRunner(const config::TimeStepperConfig& cfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::TimeStepperConfig::Type::ForwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: ForwardEuler not yet implemented");

					case config::TimeStepperConfig::Type::BackwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: BackwardEuler not yet implemented");

					case config::TimeStepperConfig::Type::GeneralizedAlpha:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: GeneralizedAlpha not yet implemented");

					case config::TimeStepperConfig::Type::RK4:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: RK4 not yet implemented");

				}

				throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: unknown timestepper type");

			}

		} // namespace timestepper
	} // namespace solver
} // namespace pdesolver

#endif
