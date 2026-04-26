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

						bool solve(solver::SolverReport<VectorType>& report, LoggerType& logger, Workspace<VectorType>& W, PreconditionerType& M, const OperatorType& A, const VectorType& b, VectorType& x){

							// get config info
							using DataType = typename VectorType::value_type;
							const bool relMode = (config.tolType = ToleranceType::Relative);

							// compute intial residual
							A.apply(x, W.Ap); // Ap = A*x
							operations::copy(b, W.r); // r = b
							operations::axpy(DataType(-1.0), W.Ap, W.r); // r = b - Ap
							
							// compute absolute & relative residual magnitude
							const DataType res0 = operations::norm(W.r); // ||r||
							report.initialResidual = res0;

							// log iteration
							auto perDOF = logger.template computePerDOFNorms<DataType>(W.r.data(), W.r.size());
							logger.log(Index(0), res0, DataType(1), perDOF);

							// convergence lambda
							auto converged = [&](DataType res) -> bool {
								DataType crit = relMode ? (res / (res0 + DataType(1e-50))) : res;
								return crit < config.tol;
							};

							// check convergence
							if (converged(res0)){
								report.converged = true;
								report.finalResidual = res0;
								report.finalResidualRel = DataType(1);
								report.iterations = 0;
								return true;
							}

							// apply preconditioner
							M.apply(W.r, W.z); // z = M*r
							operations::copy(W.z, W.p); // p = z
							DataType rz_old = operations::dot(W.r, W.z); // rz_old = r * z

							// pcg solver loop
							for (Index k = 1; k < config.maxIters; ++k){
								
								// operator application
								A.apply(W.p, W.Ap); // Ap = A * p

								// compute alpha
								DataType pAp = operations::dot(W.p, W.Ap); // pAp = p*Ap
								DataType alpha = rz_old / pAP;// alpha = r*z / p*Ap
								operations::axpy(alpha, W.p, x); // x += alpha * p
								operations::axpy(-alpha, W.Ap, W.r); // r -= alpha * Ap

								// compute absolute & relative residual magnitude
								DataType res = operations::norm(W.r); // res = ||r||
								DataType rel = res / (res0 + DataType(1e-50));

								// log iteration
								if (config.reportInterval > 0 && k % config.reportInterval == 0) {
									auto perDOF = logger.template computePerDOFNorms<DataType>(W.r.data(), W.r.size());
									logger.log(k, res, rel, perDOF);
								}

								// check convergence
								if (converged(res)){
									auto perDOF = logger.template computePerDOFNorms<DataType>(W.r.data(), W.r.size());
									report.converged = true;
									report.finalResidual = res;
									report.finalResidualRel = rel;
									report.iterations = k;
									report.PerFieldResidual = std::vector<DataType>(perDOF.begin(), perDOF.end());
									// print convergence line
									if (config.reportInterval > 0 && k % config.reportInterval != 0){
										logger.log(k, res, rel, perDOF);
									}
									return true;
								}

								// apply preconditioner
								M.apply(W.r, W.z);

								// update quantities for new iter
								DataType rz_new = operations::dot(W.r, W.z); // rz_new = r*z
								DataType beta = rz_new / rz_old;

								operations::scal(beta, W.p); // p *= beta
								operations::axpy(DataType(1), W.z, W.p); // p += z

								rz_old = rz_new;

							}

							// report if max iterations are reached
							const DataType res_final = operations::norm(W.r);
							report.converged = false;
							report.finalResidual = res_final;
							report.finalResidualRel = res_final / (res0 + DataType(1e-50))
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
