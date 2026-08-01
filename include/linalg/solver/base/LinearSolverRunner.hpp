#ifndef PDESOLVER_LINALG_SOLVER_BASE_LINEARSOLVERRUNNER_HPP
#define PDESOLVER_LINALG_SOLVER_BASE_LINEARSOLVERRUNNER_HPP

#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {

			template<typename VectorType>
			class LinearSolverRunner {
			public:

				virtual ~LinearSolverRunner() = default;

				virtual bool solve(const VectorType& b, VectorType& x, SolverReport<VectorType>& report) = 0;

			}; // class LinearSolverRunner

		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
