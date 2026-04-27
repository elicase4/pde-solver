#ifndef PDESOLVER_LINALG_CG_SOLVER_HPP
#define PDESOLVER_LINALG_CG_SOLVER_HPP

#include <cmath>
#include <cassert>
#include <vector>

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
					class Solver {
					public:
						Solver(const Config<VectorType>& cfg) : config(cfg) {}

						bool solve(solver::SolverReport<VectorType>& report, LoggerType& logger, Workspace<VectorType>& W, PreconditionerType& M, const OperatorType& A, const VectorType& b, VectorType& x);

					private:
						Config<VectorType> config;

					}; // class Solver

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#include "Solver.tpp"

#endif
