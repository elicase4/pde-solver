#ifndef PDESOLVER_LINALG_SOLVER_BASE_LINEARSOLVERALGORITHM_HPP
#define PDESOLVER_LINALG_SOLVER_BASE_LINEARSOLVERALGORITHM_HPP

#include <concepts>

#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {

			template<typename S, typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
			concept LinearSolverAlgorithm = requires(S& s, const typename S::Config& cfg, SolverReport<VectorType>& report, LoggerType& logger, typename S::Workspace& workspace, PreconditionerType& preconditioner, const OperatorType& op, const VectorType& b, VectorType& x) {
				{ S(cfg) };
				{ s.solve(report, logger, workspace, preconditioner, op, b, x) } -> std::same_as<bool>;
			}; // concept LinearSolverAlgorithm

		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
