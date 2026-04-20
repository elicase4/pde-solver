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

					template<typename OperatorType, typename VectorType, typename PreconditionerType, typename LoggerType>
					class Solver {
					public:
						Solver(const Config<VectorType>& cfg) : config(cfg) {}

						bool solve(solver::SolverReport<VectorType>& report, LoggerType& logger, Workspace<VectorType>& W, PreconditionerType& M, const OperatorType& A, const VectorType& b, VectorType& x){

							// get data type
							using DataType = typename VectorType::value_type;

							// compute intial values
							A.apply(x, W.Ap); // Ap = A*x
							operations::copy(b, W.r); // r = b
							operations::axpy(-1.0, W.Ap, W.r); // r = b - Ap
							DataType res0 = operations::norm(W.r); // ||r||
							
							// report & log intial values
							report.initialResidual = res0;
							logger.log(0, res0);

							// check convergence with initial values
							if (res0 < config.tol){
								report.converged = true;
								report.finalResidual = res0;
								report.iterations = 0;
								return true;
							}

							// apply preconditioner
							M.apply(W.r, W.z); // z = M*r
							operations::copy(W.z, W.p); // p = z

							// set previous residual product
							DataType rz_old = operations::dot(W.r, W.z); // rz_old = r * z

							// solver loop
							for (Index k = 1; k < config.maxIters; ++k){

								A.apply(W.p, W.Ap); // Ap = A * p

								// compute alpha
								DataType alpha = rz_old / operations::dot(W.p, W.Ap);// alpha = r*z / p*Ap
								operations::axpy(alpha, W.p, x); // x += alpha * p
								operations::axpy(-alpha, W.Ap, W.r); // r -= alpha * Ap

								// report & log intial values
								DataType res = operations::norm(W.r); // res = ||r||
								logger.log(k, res);

								// check convergence with initial values
								if (res < config.tol){
									report.converged = true;
									report.finalResidual = res;
									report.iterations = k;
									return true;
								}

								// apply preconditioner
								M.apply(W.r, W.z);

								// update quantities for new iter
								DataType rz_new = operations::dot(W.r, W.z); // rz_new = r*z
								DataType beta = rz_new / rz_old;

								operations::scal(beta, W.p); // p *= beta
								operations::axpy(1.0, W.z, W.p); // p += z

								rz_old = rz_new;

							}

							// set report values if max iter is reached
							report.converged = false;
							report.finalResidual = operations::norm(W.r);
							report.iterations = config.maxIters;
							
							return false;

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
