#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/TimeStepperConfig.hpp"
#include "solver/timestepper/BackwardEuler.hpp"
#include "solver/timestepper/StepSizePolicyFactory.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			template<typename StageType>
			std::unique_ptr<TimeStepperRunner> makeTimeStepperRunner(StageType& stage, const config::TimeStepperConfig& cfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::TimeStepperConfig::Type::BackwardEuler:
						return std::make_unique<BackwardEulerRunner<StageType>>(stage, cfg, makeStepSizePolicy(cfg.stepSize, equationLabel));

					case config::TimeStepperConfig::Type::ForwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: ForwardEuler not yet implemented");

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
