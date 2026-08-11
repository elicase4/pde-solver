#ifndef PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP
#define PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP

#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace solver {
		namespace nonlinear {

			// Mirrors linalg::solver::LinearSolverRunner's type-erasure role, one level up:
			// lets callers (e.g. a Stage) hold a single concrete-algorithm-agnostic pointer
			// regardless of whether NonlinearSolverConfig::Type resolved to Newton or Picard.
			//
			// No VectorType& U parameter, unlike LinearSolverRunner::solve(b, x, report) --
			// a concrete runner is constructed already bound to the NonlinearCapableStage it
			// drives, and that stage owns and updates its own solution vector internally via
			// solveLinearStep(). The runner's job is purely to orchestrate the
			// assemble/residualNorm/solveLinearStep loop against that stage.
			template<typename VectorType>
			class NonlinearSolverRunner {
			public:

				virtual ~NonlinearSolverRunner() = default;

				virtual bool solve(linalg::solver::SolverReport<VectorType>& report) = 0;

			}; // class NonlinearSolverRunner

		} // namespace nonlinear
	} // namespace solver
} // namespace pdesolver

#endif
