#ifndef PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP
#define PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP

#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace solver {
		namespace nonlinear {

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
