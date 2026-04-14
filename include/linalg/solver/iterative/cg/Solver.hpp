#ifndef PDESOLVER_LINALG_CG_SOLVER_HPP
#define PDESOLVER_LINALG_CG_SOLVER_HPP

#include <cmath>
#include <cassert>

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/solver/base/SolverReport.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorType, typename VectorType, typename PreconditionerType = void, typename LoggerType = void>
					class Solver {
					public:
						Solver(const Config<VectorType>& cfg) : config(cfg) {}

						bool solve(const OperatorType& A, const VectorType& b, VectorType& x, Workspace<VectorType>& W, SolverReport& report, PreconditionderType& M, LoggerType& logger){

							// compute intial values
							const Index n = b.size();
							A.apply(W.x, W.Ap); // Ap = A*x
							operations::copy(W.b, W.r); // b = r
							operations::axpy(-1.0, W.Ap, W.r); // r = Ap - b
							config:DataType norm0 = norm(r); // ||r||
							
							// report & log intial values
							report.initialResidual = norm0;
							logger.log(0, norm0);

							// check convergence with initial values
							if (norm0 < config.tol){
								report.converged = true;
								report.iterations = 0;
								return true;
							}

							// apply preconditioner
							M.apply(W.r, W.z); // 

						}

					private:
						Config<VectorType> config;

					}; // class Solver

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
