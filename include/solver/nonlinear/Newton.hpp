#ifndef PDESOLVER_SOLVER_NONLINEAR_NEWTON_HPP
#define PDESOLVER_SOLVER_NONLINEAR_NEWTON_HPP

#include <utility>

#include "core/Types.hpp"

#include "linalg/solver/base/SolverReport.hpp"

#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

#include "utils/logging/nonlinear/Logger.hpp"

namespace pdesolver {
	namespace solver {
		namespace nonlinear {

			template<stage::NonlinearCapableStage StageT, typename VectorType>
			class NewtonRunner : public NonlinearSolverRunner<VectorType> {
			public:

				NewtonRunner(StageT& stage, const config::NonlinearSolverConfig& cfg, utils::logging::nonlinear::Logger logger)
					: stage_(stage), cfg_(cfg), logger_(std::move(logger)) {}

				bool solve(linalg::solver::SolverReport<VectorType>& report) override {

					Real res0 = Real(0);

					for (Index iter = 0; iter < cfg_.maxIterations; ++iter) {

						stage_.assemble();
						const Real resNorm = stage_.residualNorm();

						if (iter == 0) res0 = resNorm;
						const Real resRel = (res0 > Real(0)) ? (resNorm / res0) : Real(0);

						logger_.log(iter, resNorm, resRel);

						if (resNorm < cfg_.absoluteTolerance || resRel < cfg_.relativeTolerance) {

							report.converged = true;
							report.iterations = iter;
							report.initialResidual = res0;
							report.finalResidual = resNorm;
							report.finalResidualRel = resRel;

							logger_.summary(true);
							return true;

						}

						if (!stage_.solveLinearStep()) {
							report.converged = false;
							logger_.summary(false);
							return false;
						}

					}

					report.converged = false;
					report.iterations = cfg_.maxIterations;
					logger_.summary(false);
					return false;

				}

			private:

				StageT& stage_;
				config::NonlinearSolverConfig cfg_;
				utils::logging::nonlinear::Logger logger_;

			}; // class NewtonRunner

		} // namespace nonlinear
	} // namespace solver
} // namespace pdesolver

#endif
