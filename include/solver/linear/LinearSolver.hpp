#ifndef PDESOLVER_SOLVER_LINEAR_LINEARSOLVER_HPP
#define PDESOLVER_SOLVER_LINEAR_LINEARSOLVER_HPP

#include <concepts>

namespace pdesolver {
	namespace solver {

		template<typename S, typename OperatorType, typename VectorType, typename ReportType>
		concept LinearSolver = requires(S& s, const OperatorType& op, const VectorType& F, VectorType& U, ReportType report) {
			{ s.solve(op, F, U, report) } -> std::same_as<bool>;
		}; // concept LinearSolver

	} // namespace solver
} // namespace pdesolver

#endif
