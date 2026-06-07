#ifndef PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVER_HPP
#define PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVER_HPP

#include <concepts>

namespace pdesolver {
	namespace solver {

		template<typename S, typename TangentOperatorType, typename ResidualOperatorType, typename VectorType, typename ReportType>
		concept NonlinearSolver = requires(S& s, const TangentOperatorType& J, const ResidualOperatorType& R, VectorType& U, ReportType report) {
			{ s.solve(J, R, U, report) } -> std::same_as<bool>;
		}; // concept NonlinearSolver

	} // namespace solver
} // namespace pdesolver

#endif
