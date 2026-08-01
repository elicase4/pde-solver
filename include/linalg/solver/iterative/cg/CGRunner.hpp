#ifndef PDESOLVER_LINALG_SOLVER_ITERATIVE_CG_CGRUNNER_HPP
#define PDESOLVER_LINALG_SOLVER_ITERATIVE_CG_CGRUNNER_HPP

#include "core/Types.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/base/LinearSolverAlgorithm.hpp"

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
					class CGRunner : public linalg::solver::LinearSolverRunner<VectorType> {
					public:

						using AlgorithmT = Solver<OperatorType, VectorType, PreconditionerType, LoggerType>;

						static_assert( linalg::solver::LinearSolverAlgorithm<AlgorithmT, OperatorType, VectorType, PreconditionerType, LoggerType>, "CGRunner: cg::Solver no longer satisfies linalg::solver::LinearSolverAlgorithm. Check its Config/Workspace aliases and solve() signature");

						CGRunner(const OperatorType& op, Index n, const typename AlgorithmT::Config& cfg, LoggerType logger);

						bool solve(const VectorType& b, VectorType& x, linalg::solver::SolverReport<VectorType>& report) override;

					private:

						OperatorType op_;
						typename AlgorithmT::Workspace workspace_;
						PreconditionerType preconditioner_;
						LoggerType logger_;
						AlgorithmT solver_;

					}; // class CGRunner

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#include "linalg/solver/iterative/cg/CGRunner.tpp"

#endif
