#ifndef PDESOLVER_SOLVER_TIMESTEPPER_BACKWARDEULER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_BACKWARDEULER_HPP

#include <algorithm>
#include <memory>

#include "core/Types.hpp"

#include "solver/config/TimeStepperConfig.hpp"
#include "solver/stage/TransientCapableStage.hpp"
#include "solver/timestepper/TimeStepSizePolicy.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

#include "utils/logging/timestepper/Logger.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			template<stage::TransientCapableStage StageType>
			class BackwardEulerRunner : public TimeStepperRunner {
			public:

				BackwardEulerRunner(StageType& stage, const config::TimeStepperConfig& cfg, std::unique_ptr<TimeStepSizePolicy> policy, utils::logging::timestepper::Logger logger)
					: stage_(stage), policy_(std::move(policy)), logger_(std::move(logger)), time_(cfg.t0), tf_(cfg.tf) {}

				// summary is printed here, once this runner's work is over, rather than by the
				// outer Transient driver -- mirrors CGRunner calling logger.summary(...) itself
				// at the point its own work concludes.
				~BackwardEulerRunner() override { logger_.summary(finished()); }

				bool step() override {

					Real dt;
					Index attempts = 0;

					do {

						++attempts;

						// never overshoot tf
						dt = std::min(policy_->dt(), tf_ - time_);

						stage_.setDt(dt);
						stage_.setTime(time_ + dt);
						stage_.assemble();

						if (!stage_.solve()) {
							return false;
						}

						policy_->onStepComplete(stage_.residualNorm());

					} while (policy_->rejectLastStep());

					stage_.advance();
					time_ += dt;
					++step_;
					stage_.onStepComplete(step_, time_);

					logger_.log(step_, time_, dt, attempts, stage_.residualNorm());

					return true;

				}

				bool finished() const override { return time_ >= tf_ - Real(1e-9); }

				Index currentStep() const override { return step_; }

				Real currentTime() const override { return time_; }

			private:

				StageType& stage_;
				std::unique_ptr<TimeStepSizePolicy> policy_;
				utils::logging::timestepper::Logger logger_;

				Real time_;
				Real tf_;
				Index step_ = 0;

			}; // class BackwardEulerRunner

		} // namespace timestepper
	} // namespace solver
} // namespace pdesolver

#endif
