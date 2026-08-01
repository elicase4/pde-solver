#ifndef PDESOLVER_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP
#define PDESOLVER_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/Types.hpp"

#include "solver/config/LinearSolverConfig.hpp"

#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/iterative/cg/CGRunner.hpp"
#include "linalg/solver/iterative/cg/Config.hpp"

#include "linalg/solver/preconditioner/Identity.hpp"

#include "utils/logging/solver/ConsoleLogger.hpp"

namespace pdesolver {
	namespace solver {
		namespace linear {

			template<typename OperatorType, typename VectorType>
			std::unique_ptr<linalg::solver::LinearSolverRunner<VectorType>> makeLinearSolverRunner(const OperatorType& op, Index n, const config::LinearSolverConfig& cfg, const std::string& equationLabel) {

				if (cfg.preconditioner.type != config::PreconditionerConfig::Type::Identity) {
					throw std::runtime_error("LinearSolverFactory: only the identity preconditioner is implemented so far");
				}

				using PreconditionerT = linalg::solver::preconditioner::Identity<VectorType>;
				using LoggerT = utils::logging::ConsoleLogger;

				switch (cfg.type) {

					case config::LinearSolverConfig::Type::CG: {

						using CGConfigT = linalg::solver::iterative::cg::Config<VectorType>;
						CGConfigT cgCfg;
						cgCfg.tol = static_cast<typename CGConfigT::DataType>(cfg.tolerance);
						cgCfg.maxIters = cfg.maxIterations;

						LoggerT logger("CG", equationLabel);

						return std::make_unique<linalg::solver::iterative::cg::CGRunner<OperatorType, VectorType, PreconditionerT, LoggerT>>(op, n, cgCfg, std::move(logger));

					}

					case config::LinearSolverConfig::Type::GMRES:
						throw std::runtime_error("LinearSolverFactory: GMRES not yet implemented");

					case config::LinearSolverConfig::Type::BiCGSTAB:
						throw std::runtime_error("LinearSolverFactory: BiCGSTAB not yet implemented");

					case config::LinearSolverConfig::Type::LU:
						throw std::runtime_error("LinearSolverFactory: LU not yet implemented");

				}

				throw std::runtime_error("LinearSolverFactory: unknown linear solver type");

			}

		} // namespace linear
	} // namespace solver
} // namespace pdesolver

#endif
