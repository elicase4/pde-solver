#ifndef PDESOLVER_APPLICATION_HEATEQ_STAGE_HEATSTAGE_HPP
#define PDESOLVER_APPLICATION_HEATEQ_STAGE_HEATSTAGE_HPP

#include <memory>

#include "core/Types.hpp"

#include "linalg/solver/base/SolverReport.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/nonlinear/NonlinearSolverFactory.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace stage {

				template<typename ProblemT>
				class HeatStage {
				public:

					explicit HeatStage(ProblemT& problem) : problem_(problem) {

						if (solver::isNonlinear(problem_.solverInstance().mode)) {
							// VectorType isn't deducible from any parameter here (only StageT is,
							// from the stage argument) -- both must be given explicitly.
							nonlinearSolverRunner_ = solver::nonlinear::makeNonlinearSolverRunner<HeatStage, typename ProblemT::VectorT>(
								*this, *problem_.solverInstance().nonlinear, "heateq");
						}

					}

					void initialize() {}

					void assemble() { problem_.assembleSystem(currentTime_); }

					bool solve() {

						if (nonlinearSolverRunner_) {
							linalg::solver::SolverReport<typename ProblemT::VectorT> report;
							return nonlinearSolverRunner_->solve(report);
						}

						return problem_.solveLinear();

					}

					void finalize() {}

					// NonlinearCapableStage forwarding -- only ever called by
					// nonlinearSolverRunner_'s (still unimplemented) internal loop.
					Real residualNorm() const { return problem_.residualNorm(); }

					bool solveLinearStep() { return problem_.solveLinearStep(); }

					decltype(auto) solution() const { return problem_.solution(); }

				private:

					ProblemT& problem_;

					Real currentTime_ = 0.0;

					std::unique_ptr<solver::nonlinear::NonlinearSolverRunner<typename ProblemT::VectorT>> nonlinearSolverRunner_;

				}; // class HeatStage

			} // namespace stage
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
