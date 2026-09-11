#ifndef PDESOLVER_SOLVER_TIMESTEPPER_BACKWARDEULER_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_BACKWARDEULER_HPP

#include <algorithm>
#include <memory>

#include "core/Types.hpp"

#include "solver/config/TimeStepperConfig.hpp"
#include "solver/stage/TransientCapableStage.hpp"
#include "solver/timestepper/TimeStepSizePolicy.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			template<stage::TransientCapableStage StageType>
			class BackwardEulerRunner : public TimeStepperRunner {
			public:

				BackwardEulerRunner(StageType& stage, const config::TimeStepperConfig& cfg, std::unique_ptr<TimeStepSizePolicy> policy)
					: stage_(stage), policy_(std::move(policy)), time_(cfg.t0), tf_(cfg.tf) {}

				bool step() override {

					Real dt;

					do {

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

					return true;

				}

				bool finished() const override { return time_ >= tf_ - Real(1e-9); }

				Index currentStep() const override { return step_; }

				Real currentTime() const override { return time_; }

			private:

				StageType& stage_;
				std::unique_ptr<TimeStepSizePolicy> policy_;

				Real time_;
				Real tf_;
				Index step_ = 0;

			}; // class BackwardEulerRunner

		} // namespace timestepper
	} // namespace solver
} // namespace pdesolver

#endif
